/**
 * @fileoverview Client-side BAM header extraction.
 *
 * Reads just enough of a BAM file to extract the SAM header text,
 * avoiding the need to upload the entire file.
 * @module utils/headerExtractor
 */

import { MAX_BAM_HEADER_READ } from './helpers.js';

/**
 * Concatenate an array of Uint8Array chunks into a single Uint8Array.
 * @param {Uint8Array[]} chunks
 * @returns {Uint8Array}
 */
function concatChunks(chunks) {
    const totalLength = chunks.reduce((sum, c) => sum + c.length, 0);
    const result = new Uint8Array(totalLength);
    let offset = 0;
    for (const chunk of chunks) {
        result.set(chunk, offset);
        offset += chunk.length;
    }
    return result;
}

/**
 * Read a slice of a File as an ArrayBuffer.
 * @param {File} file
 * @param {number} start
 * @param {number} end
 * @returns {Promise<ArrayBuffer>}
 */
function readFileSlice(file, start, end) {
    return new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = () => resolve(reader.result);
        reader.onerror = () => reject(new Error('Failed to read file'));
        reader.readAsArrayBuffer(file.slice(start, end));
    });
}

/**
 * Decompress a single gzip member (BGZF block) using DecompressionStream.
 * @param {Uint8Array} block - Complete gzip member
 * @returns {Promise<Uint8Array>} Decompressed data
 */
async function decompressGzipBlock(block) {
    const ds = new DecompressionStream('gzip');
    const writer = ds.writable.getWriter();
    const reader = ds.readable.getReader();

    try {
        await writer.write(block);
        await writer.close();

        const chunks = [];
        while (true) {
            const { done, value } = await reader.read();
            if (done) break;
            chunks.push(value);
        }

        return concatChunks(chunks);
    } catch (err) {
        reader.cancel().catch(() => {});
        writer.abort(err).catch(() => {});
        throw err;
    }
}

/**
 * Decompress consecutive BGZF blocks, stopping as soon as we have
 * enough uncompressed data for the BAM header.
 * @param {Uint8Array} bytes - Raw BGZF data
 * @param {number} [neededBytes=0] - Stop after accumulating this many uncompressed bytes (0 = all)
 * @returns {Promise<Uint8Array>} Concatenated decompressed data
 */
async function decompressBgzfBlocks(bytes, neededBytes = 0) {
    const chunks = [];
    let offset = 0;
    let accumulated = 0;

    while (offset < bytes.length) {
        // Check for gzip magic
        if (offset + 18 > bytes.length || bytes[offset] !== 0x1f || bytes[offset + 1] !== 0x8b) {
            break;
        }

        // FLG must have FEXTRA (0x04)
        if ((bytes[offset + 3] & 0x04) === 0) break;

        const view = new DataView(bytes.buffer, bytes.byteOffset + offset, bytes.byteLength - offset);
        const xlen = view.getUint16(10, true);
        const extraEnd = 12 + xlen;

        // Find BC subfield to get BSIZE
        let bsize = -1;
        let pos = 12;
        while (pos + 4 <= extraEnd) {
            if (bytes[offset + pos] === 0x42 && bytes[offset + pos + 1] === 0x43) {
                const slen = view.getUint16(pos + 2, true);
                if (slen === 2) {
                    bsize = view.getUint16(pos + 4, true) + 1;
                    break;
                }
            }
            const slen = view.getUint16(pos + 2, true);
            pos += 4 + slen;
        }

        if (bsize <= 0 || offset + bsize > bytes.length) break;

        const block = bytes.subarray(offset, offset + bsize);
        try {
            const decompressed = await decompressGzipBlock(block);
            if (decompressed.length === 0) break; // EOF block
            chunks.push(decompressed);
            accumulated += decompressed.length;
        } catch {
            break;
        }

        offset += bsize;

        // Stop early once we have enough data for the header
        if (neededBytes > 0 && accumulated >= neededBytes) break;
    }

    return concatChunks(chunks);
}

/**
 * Extract the SAM header text from a BAM file.
 * @param {File} file - BAM file
 * @returns {Promise<{header: string, filename: string}>}
 * @throws {Error} If the file is not a valid BAM or extraction fails
 */
export async function extractBamHeader(file) {
    const readSize = Math.min(file.size, MAX_BAM_HEADER_READ);
    const buffer = await readFileSlice(file, 0, readSize);
    const bytes = new Uint8Array(buffer);

    // First pass: decompress enough to read the l_text field (8 bytes minimum)
    let uncompressed = await decompressBgzfBlocks(bytes, 8);

    // Verify BAM magic: "BAM\1"
    if (uncompressed.length < 8 ||
        uncompressed[0] !== 0x42 || uncompressed[1] !== 0x41 ||
        uncompressed[2] !== 0x4d || uncompressed[3] !== 0x01) {
        throw new Error('Not a valid BAM file (bad magic bytes)');
    }

    const view = new DataView(uncompressed.buffer);
    const headerLength = view.getInt32(4, true);
    if (headerLength < 0) {
        throw new Error('BAM header length is negative');
    }

    // If we don't have enough data yet, decompress more blocks
    const needed = 8 + headerLength;
    if (uncompressed.length < needed) {
        uncompressed = await decompressBgzfBlocks(bytes, needed);
        if (uncompressed.length < needed) {
            throw new Error('BAM header length exceeds available data');
        }
    }

    const decoder = new TextDecoder('ascii');
    const headerText = decoder.decode(uncompressed.subarray(8, 8 + headerLength));

    return { header: headerText.trimEnd(), filename: file.name };
}

/**
 * Check if a filename indicates a BAM file.
 * @param {string} filename
 * @returns {boolean}
 */
export function isBamFile(filename) {
    return filename.toLowerCase().endsWith('.bam');
}
