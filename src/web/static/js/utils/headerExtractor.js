/**
 * @fileoverview Client-side BAM header extraction.
 *
 * Reads just enough of a BAM file to extract the SAM header text,
 * avoiding the need to upload the entire file.
 * @module utils/headerExtractor
 */

/** Max bytes to read from file for header extraction. */
const MAX_HEADER_READ = 4 * 1024 * 1024;

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

        const totalLength = chunks.reduce((sum, c) => sum + c.length, 0);
        const result = new Uint8Array(totalLength);
        let offset = 0;
        for (const chunk of chunks) {
            result.set(chunk, offset);
            offset += chunk.length;
        }
        return result;
    } catch (err) {
        reader.cancel().catch(() => {});
        writer.abort(err).catch(() => {});
        throw err;
    }
}

/**
 * Decompress consecutive BGZF blocks from raw bytes.
 * @param {Uint8Array} bytes - Raw BGZF data
 * @returns {Promise<Uint8Array>} Concatenated decompressed data
 */
async function decompressBgzfBlocks(bytes) {
    const chunks = [];
    let offset = 0;

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
        } catch {
            break;
        }

        offset += bsize;
    }

    const totalLength = chunks.reduce((sum, c) => sum + c.length, 0);
    const result = new Uint8Array(totalLength);
    let writeOffset = 0;
    for (const chunk of chunks) {
        result.set(chunk, writeOffset);
        writeOffset += chunk.length;
    }
    return result;
}

/**
 * Extract the SAM header text from a BAM file.
 * @param {File} file - BAM file
 * @returns {Promise<{header: string, filename: string}>}
 * @throws {Error} If the file is not a valid BAM or extraction fails
 */
export async function extractBamHeader(file) {
    const readSize = Math.min(file.size, MAX_HEADER_READ);
    const buffer = await readFileSlice(file, 0, readSize);
    const bytes = new Uint8Array(buffer);

    const uncompressed = await decompressBgzfBlocks(bytes);
    const view = new DataView(uncompressed.buffer);

    // Verify BAM magic: "BAM\1"
    if (uncompressed.length < 8 ||
        uncompressed[0] !== 0x42 || uncompressed[1] !== 0x41 ||
        uncompressed[2] !== 0x4d || uncompressed[3] !== 0x01) {
        throw new Error('Not a valid BAM file (bad magic bytes)');
    }

    const headerLength = view.getInt32(4, true);
    if (headerLength < 0 || headerLength > uncompressed.length - 8) {
        throw new Error('BAM header length exceeds available data');
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

/**
 * Check if a filename indicates a CRAM file.
 * @param {string} filename
 * @returns {boolean}
 */
export function isCramFile(filename) {
    return filename.toLowerCase().endsWith('.cram');
}
