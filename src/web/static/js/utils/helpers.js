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
@SQ\tSN:chr2\tLN:242193529\tM5:f98db672eb0993dcfdabafe2a882905c
@SQ\tSN:chr3\tLN:198295559\tM5:76635a41ea913a405ded820447d067b0
@SQ\tSN:chr4\tLN:190214555\tM5:3210fecf1eb92d5489da4346b3fddc6e
@SQ\tSN:chr5\tLN:181538259\tM5:a811b3dc9fe66af729dc0dddf7fa4f13
@SQ\tSN:chr6\tLN:170805979\tM5:5691468a67c7e7a7b5f2a3a683792c29
@SQ\tSN:chr7\tLN:159345973\tM5:cc044cc2256a1141212660fb07b6171e
@SQ\tSN:chr8\tLN:145138636\tM5:c67955b5f7815a9a1edfaa15893d3616
@SQ\tSN:chr9\tLN:138394717\tM5:6c198acf68b5af7b9d676dfdd531b5de
@SQ\tSN:chr10\tLN:133797422\tM5:c0eeee7acfdaf31b770a509bdaa6e51a
@SQ\tSN:chr11\tLN:135086622\tM5:1511375dc2dd1b633af8cf439ae90cec
@SQ\tSN:chr12\tLN:133275309\tM5:96e414eace405d8c27a6d35ba19df56f
@SQ\tSN:chr13\tLN:114364328\tM5:a5437debe2ef9c9ef8f3ea2874ae1d82
@SQ\tSN:chr14\tLN:107043718\tM5:e0f0eecc3bcab6178c62b6211565c807
@SQ\tSN:chr15\tLN:101991189\tM5:f036bd11158407596ca6bf3581454706
@SQ\tSN:chr16\tLN:90338345\tM5:db2d37c8b7d019caaf2dd64ba3a6f33a
@SQ\tSN:chr17\tLN:83257441\tM5:f9a0fb01553adb183568e3eb9d8626db
@SQ\tSN:chr18\tLN:80373285\tM5:11eeaa801f6b0e2e36a1138616b8ee9a
@SQ\tSN:chr19\tLN:58617616\tM5:85f9f4fc152c58cb7913c06d6b98573a
@SQ\tSN:chr20\tLN:64444167\tM5:b18e6c531b0bd70e949a7fc20859cb01
@SQ\tSN:chr21\tLN:46709983\tM5:974dc7aec0b755b19f031418fdedf293
@SQ\tSN:chr22\tLN:50818468\tM5:ac37ec46683600f808cdd41eac1d55cd
@SQ\tSN:chrX\tLN:156040895\tM5:2b3a55ff7f58eb308420c8a9b11cac50
@SQ\tSN:chrY\tLN:57227415\tM5:ce3e31103314a704255f3cd90369ecce
@SQ\tSN:chrM\tLN:16569\tM5:c68f52674c9fb33aef52dcf399755519`,
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
4\tassembled-molecule\t4\tChromosome\tCM000666.2\t=\tNC_000004.12\tPrimary Assembly\t190214555\tchr4
5\tassembled-molecule\t5\tChromosome\tCM000667.2\t=\tNC_000005.10\tPrimary Assembly\t181538259\tchr5
6\tassembled-molecule\t6\tChromosome\tCM000668.2\t=\tNC_000006.12\tPrimary Assembly\t170805979\tchr6
7\tassembled-molecule\t7\tChromosome\tCM000669.2\t=\tNC_000007.14\tPrimary Assembly\t159345973\tchr7
8\tassembled-molecule\t8\tChromosome\tCM000670.2\t=\tNC_000008.11\tPrimary Assembly\t145138636\tchr8
9\tassembled-molecule\t9\tChromosome\tCM000671.2\t=\tNC_000009.12\tPrimary Assembly\t138394717\tchr9
10\tassembled-molecule\t10\tChromosome\tCM000672.2\t=\tNC_000010.11\tPrimary Assembly\t133797422\tchr10
11\tassembled-molecule\t11\tChromosome\tCM000673.2\t=\tNC_000011.10\tPrimary Assembly\t135086622\tchr11
12\tassembled-molecule\t12\tChromosome\tCM000674.2\t=\tNC_000012.12\tPrimary Assembly\t133275309\tchr12
13\tassembled-molecule\t13\tChromosome\tCM000675.2\t=\tNC_000013.11\tPrimary Assembly\t114364328\tchr13
14\tassembled-molecule\t14\tChromosome\tCM000676.2\t=\tNC_000014.9\tPrimary Assembly\t107043718\tchr14
15\tassembled-molecule\t15\tChromosome\tCM000677.2\t=\tNC_000015.10\tPrimary Assembly\t101991189\tchr15
16\tassembled-molecule\t16\tChromosome\tCM000678.2\t=\tNC_000016.10\tPrimary Assembly\t90338345\tchr16
17\tassembled-molecule\t17\tChromosome\tCM000679.2\t=\tNC_000017.11\tPrimary Assembly\t83257441\tchr17
18\tassembled-molecule\t18\tChromosome\tCM000680.2\t=\tNC_000018.10\tPrimary Assembly\t80373285\tchr18
19\tassembled-molecule\t19\tChromosome\tCM000681.2\t=\tNC_000019.10\tPrimary Assembly\t58617616\tchr19
20\tassembled-molecule\t20\tChromosome\tCM000682.2\t=\tNC_000020.11\tPrimary Assembly\t64444167\tchr20
21\tassembled-molecule\t21\tChromosome\tCM000683.2\t=\tNC_000021.9\tPrimary Assembly\t46709983\tchr21
22\tassembled-molecule\t22\tChromosome\tCM000684.2\t=\tNC_000022.11\tPrimary Assembly\t50818468\tchr22
X\tassembled-molecule\tX\tChromosome\tCM000685.2\t=\tNC_000023.11\tPrimary Assembly\t156040895\tchrX
Y\tassembled-molecule\tY\tChromosome\tCM000686.2\t=\tNC_000024.10\tPrimary Assembly\t57227415\tchrY
MT\tassembled-molecule\tMT\tMitochondrion\tJ01415.2\t=\tNC_012920.1\tnon-nuclear\t16569\tchrM`,
    'vcf': `##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422,assembly=GRCh38,md5=6aef897c3d6ff0c78aff06ac189178dd>
##contig=<ID=chr2,length=242193529,assembly=GRCh38,md5=f98db672eb0993dcfdabafe2a882905c>
##contig=<ID=chr3,length=198295559,assembly=GRCh38,md5=76635a41ea913a405ded820447d067b0>
##contig=<ID=chr4,length=190214555,assembly=GRCh38,md5=3210fecf1eb92d5489da4346b3fddc6e>
##contig=<ID=chr5,length=181538259,assembly=GRCh38,md5=a811b3dc9fe66af729dc0dddf7fa4f13>
##contig=<ID=chr6,length=170805979,assembly=GRCh38,md5=5691468a67c7e7a7b5f2a3a683792c29>
##contig=<ID=chr7,length=159345973,assembly=GRCh38,md5=cc044cc2256a1141212660fb07b6171e>
##contig=<ID=chr8,length=145138636,assembly=GRCh38,md5=c67955b5f7815a9a1edfaa15893d3616>
##contig=<ID=chr9,length=138394717,assembly=GRCh38,md5=6c198acf68b5af7b9d676dfdd531b5de>
##contig=<ID=chr10,length=133797422,assembly=GRCh38,md5=c0eeee7acfdaf31b770a509bdaa6e51a>
##contig=<ID=chr11,length=135086622,assembly=GRCh38,md5=1511375dc2dd1b633af8cf439ae90cec>
##contig=<ID=chr12,length=133275309,assembly=GRCh38,md5=96e414eace405d8c27a6d35ba19df56f>
##contig=<ID=chr13,length=114364328,assembly=GRCh38,md5=a5437debe2ef9c9ef8f3ea2874ae1d82>
##contig=<ID=chr14,length=107043718,assembly=GRCh38,md5=e0f0eecc3bcab6178c62b6211565c807>
##contig=<ID=chr15,length=101991189,assembly=GRCh38,md5=f036bd11158407596ca6bf3581454706>
##contig=<ID=chr16,length=90338345,assembly=GRCh38,md5=db2d37c8b7d019caaf2dd64ba3a6f33a>
##contig=<ID=chr17,length=83257441,assembly=GRCh38,md5=f9a0fb01553adb183568e3eb9d8626db>
##contig=<ID=chr18,length=80373285,assembly=GRCh38,md5=11eeaa801f6b0e2e36a1138616b8ee9a>
##contig=<ID=chr19,length=58617616,assembly=GRCh38,md5=85f9f4fc152c58cb7913c06d6b98573a>
##contig=<ID=chr20,length=64444167,assembly=GRCh38,md5=b18e6c531b0bd70e949a7fc20859cb01>
##contig=<ID=chr21,length=46709983,assembly=GRCh38,md5=974dc7aec0b755b19f031418fdedf293>
##contig=<ID=chr22,length=50818468,assembly=GRCh38,md5=ac37ec46683600f808cdd41eac1d55cd>
##contig=<ID=chrX,length=156040895,assembly=GRCh38,md5=2b3a55ff7f58eb308420c8a9b11cac50>
##contig=<ID=chrY,length=57227415,assembly=GRCh38,md5=ce3e31103314a704255f3cd90369ecce>
##contig=<ID=chrM,length=16569,assembly=GRCh38,md5=c68f52674c9fb33aef52dcf399755519>`
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