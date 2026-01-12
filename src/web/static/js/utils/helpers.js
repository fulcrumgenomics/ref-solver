/**
 * @fileoverview Utility helpers and shared constants
 * @module utils/helpers
 */

/**
 * Example data for demonstrations
 * @type {Object.<string, string>}
 */
export const examples = {
    'sam-dict': `@SQ\tSN:chr1\tLN:248956422\tM5:6aef897c3d6ff0c78aff06ac189178dd
@SQ\tSN:chr2\tLN:242193529\tM5:f98db672eb0993dcfdabafe2a882905c`,
    'assembly-report': `# Assembly name:  GRCh38.p14
# Organism name:  Homo sapiens (human)
# Infraspecific name:  strain=reference
# Taxid:          9606
# BioSample:      SAMN02981242
# BioProject:     PRJNA31257
# Submitter:      Genome Reference Consortium
# Date:           2013-12-17
# Assembly type:  haploid-with-alt-loci
# Release type:   patch
# Assembly level: Chromosome
# Genome representation: full
#
# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name
1\tassembled-molecule\t1\tChromosome\tCM000663.2\t=\tNC_000001.11\tPrimary Assembly\t248956422\tchr1
2\tassembled-molecule\t2\tChromosome\tCM000664.2\t=\tNC_000002.12\tPrimary Assembly\t242193529\tchr2
3\tassembled-molecule\t3\tChromosome\tCM000665.2\t=\tNC_000003.12\tPrimary Assembly\t198295559\tchr3
X\tassembled-molecule\tX\tChromosome\tCM000685.2\t=\tNC_000023.11\tPrimary Assembly\t156040895\tchrX
Y\tassembled-molecule\tY\tChromosome\tCM000686.2\t=\tNC_000024.10\tPrimary Assembly\t57227415\tchrY`,
    'vcf': `##contig=<ID=chr1,length=248956422,assembly=b38,md5=6aef897c3d6ff0c78aff06ac189178dd,species="Homo sapiens",taxonomy=x>
##contig=<ID=chr2,length=242193529,assembly=b38,md5=f98db672eb0993dcfdabafe2a882905c,species="Homo sapiens",taxonomy=x>
##contig=<ID=chr3,length=198295559,assembly=b38,md5=76635a41ea913a405ded8c937d83cb04,species="Homo sapiens",taxonomy=x>
##contig=<ID=chr4,length=190214555,assembly=b38,md5=3210fecf1eb92d5489da4346b3fddc6e,species="Homo sapiens",taxonomy=x>
##contig=<ID=chr5,length=181538259,assembly=b38,md5=a811b3dc9fe66af729dc0dddf7fa4f33,species="Homo sapiens",taxonomy=x>
##contig=<ID=chrX,length=156040895,assembly=b38,md5=2b3a55ff7f58eb308420c8a9df3f0b28,species="Homo sapiens",taxonomy=x>
##contig=<ID=chrY,length=57227415,assembly=b38,md5=ce3e31103314a704255595d80f2b6d32,species="Homo sapiens",taxonomy=x>`
};

/**
 * Maximum file size for uploads (16MB)
 * @constant {number}
 */
export const MAX_FILE_SIZE = 16 * 1024 * 1024;

/**
 * Maximum file size for text files (1MB)
 * @constant {number}
 */
export const MAX_TEXT_FILE_SIZE = 1024 * 1024;

/**
 * Debounce delay in milliseconds
 * @constant {number}
 */
export const DEBOUNCE_DELAY = 300;

/**
 * Escape HTML to prevent XSS attacks
 * @param {string|null|undefined} text - Text to escape
 * @returns {string} Escaped HTML string
 */
export function escapeHtml(text) {
    if (!text) return '';
    const div = document.createElement('div');
    div.textContent = text;
    return div.innerHTML;
}

/**
 * Validate that a value is within a numeric range
 * @param {number} value - Value to validate
 * @param {number} min - Minimum allowed value
 * @param {number} max - Maximum allowed value
 * @returns {number} Clamped value within range
 */
export function clamp(value, min, max) {
    return Math.min(Math.max(value, min), max);
}

/**
 * Validate file size against maximum allowed
 * @param {File} file - File to validate
 * @param {number} maxSize - Maximum allowed size in bytes
 * @returns {{valid: boolean, error?: string}} Validation result
 */
export function validateFileSize(file, maxSize) {
    if (!file) {
        return { valid: false, error: 'No file provided' };
    }
    if (file.size > maxSize) {
        const maxMB = (maxSize / (1024 * 1024)).toFixed(1);
        const fileMB = (file.size / (1024 * 1024)).toFixed(1);
        return {
            valid: false,
            error: `File size (${fileMB}MB) exceeds maximum allowed (${maxMB}MB)`
        };
    }
    return { valid: true };
}

/**
 * Create a debounced version of a function
 * @template {Function} T
 * @param {T} func - Function to debounce
 * @param {number} wait - Debounce delay in milliseconds
 * @returns {T} Debounced function
 */
export function debounce(func, wait) {
    let timeout;
    return function executedFunction(...args) {
        const later = () => {
            clearTimeout(timeout);
            func(...args);
        };
        clearTimeout(timeout);
        timeout = setTimeout(later, wait);
    };
}

/**
 * Format file size for display
 * @param {number} bytes - Size in bytes
 * @returns {string} Formatted size string
 */
export function formatFileSize(bytes) {
    if (bytes === 0) return '0 Bytes';
    const k = 1024;
    const sizes = ['Bytes', 'KB', 'MB', 'GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
}