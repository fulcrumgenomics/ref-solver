/**
 * @fileoverview Tab Manager for input format switching
 * @module managers/TabManager
 */

/**
 * @typedef {'sam-dict'|'assembly-report'|'vcf'|'binary'} TabId
 * Valid tab identifiers
 */

/**
 * @typedef {'text'|'assembly'|'vcf'|'binary'} FileFormat
 * Supported file format types
 */

/**
 * @typedef {Object} TextInput
 * @property {'text'} type - Input type
 * @property {string} content - Text content
 * @property {string} format - Format identifier (sam, assembly, vcf)
 * @property {string} [filename] - Optional filename
 */

/**
 * @typedef {Object} FileInput
 * @property {'file'} type - Input type
 * @property {File} file - File object
 * @property {FileFormat} format - File format
 */

/**
 * @typedef {TextInput|FileInput|null} CurrentInput
 * Current input data for form submission
 */

/**
 * Tab Manager
 * Handles tab switching and input validation for different file formats
 * @class
 */
export class TabManager {
    /**
     * Creates a new TabManager instance
     * @constructor
     */
    constructor() {
        /** @type {TabId} Currently active tab */
        this.activeTab = 'sam-dict';

        /** @type {File|null} Currently selected file */
        this.currentFile = null;

        /** @type {FileFormat} Current file format */
        this.currentFormat = 'text';
    }

    /**
     * Clear the current file selection
     * @returns {void}
     */
    clearCurrentFile() {
        this.currentFile = null;
        this.currentFormat = 'text';
        // Clear any file previews
        const preview = document.getElementById('binary-preview');
        if (preview) {
            preview.style.display = 'none';
        }
    }

    /**
     * Switch to a different tab
     * @param {TabId} tabId - ID of the tab to switch to
     * @returns {void}
     */
    switchTab(tabId) {
        // Validate tab ID
        const validTabs = ['sam-dict', 'assembly-report', 'vcf', 'binary'];
        if (!validTabs.includes(tabId)) {
            console.warn(`Invalid tab ID: ${tabId}`);
            return;
        }

        // Hide all tab contents
        document.querySelectorAll('.tab-content').forEach(tab => {
            tab.classList.remove('active');
        });

        // Remove active class from all tab buttons
        document.querySelectorAll('.tab-button').forEach(btn => {
            btn.classList.remove('active');
        });

        // Show selected tab content
        const tabContent = document.getElementById('tab-' + tabId);
        if (tabContent) {
            tabContent.classList.add('active');
        }

        // Find and activate the correct tab button
        document.querySelectorAll('.tab-button').forEach(btn => {
            if (btn.onclick && btn.onclick.toString().includes(`'${tabId}'`)) {
                btn.classList.add('active');
            }
        });

        this.activeTab = tabId;
        this.clearCurrentFile();
        this.validateFormat();
    }

    /**
     * Validate current input and update submit button state
     * @returns {void}
     */
    validateFormat() {
        const submitBtn = document.getElementById('submit-btn');
        if (submitBtn) {
            submitBtn.disabled = !this.hasValidContent();
        }
    }

    /**
     * Check if the current tab has valid content
     * @returns {boolean} True if valid content exists
     */
    hasValidContent() {
        switch (this.activeTab) {
            case 'sam-dict': {
                const input = document.getElementById('header-input');
                return (input && input.value.trim() !== '') || this.currentFile !== null;
            }
            case 'assembly-report': {
                const input = document.getElementById('assembly-input');
                return (input && input.value.trim() !== '') || this.currentFile !== null;
            }
            case 'vcf': {
                const input = document.getElementById('vcf-input');
                return (input && input.value.trim() !== '') || this.currentFile !== null;
            }
            case 'binary':
                return this.currentFile !== null;
            default:
                return false;
        }
    }

    /**
     * Get the current input data for form submission
     * @returns {CurrentInput} Input data object or null if no valid input
     */
    getCurrentInput() {
        switch (this.activeTab) {
            case 'sam-dict': {
                const textarea = document.getElementById('header-input');
                const textContent = textarea ? textarea.value.trim() : '';
                if (textContent) {
                    return {
                        type: 'text',
                        content: textContent,
                        format: 'sam'
                    };
                } else if (this.currentFile) {
                    return {
                        type: 'file',
                        file: this.currentFile,
                        format: this.currentFormat
                    };
                }
                break;
            }
            case 'assembly-report': {
                const textarea = document.getElementById('assembly-input');
                const textContent = textarea ? textarea.value.trim() : '';
                if (textContent) {
                    return {
                        type: 'text',
                        content: textContent,
                        format: 'assembly'
                    };
                } else if (this.currentFile) {
                    return {
                        type: 'file',
                        file: this.currentFile,
                        format: this.currentFormat
                    };
                }
                break;
            }
            case 'vcf': {
                const textarea = document.getElementById('vcf-input');
                const textContent = textarea ? textarea.value.trim() : '';
                if (textContent) {
                    return {
                        type: 'text',
                        content: textContent,
                        format: 'vcf'
                    };
                } else if (this.currentFile) {
                    return {
                        type: 'file',
                        file: this.currentFile,
                        format: this.currentFormat
                    };
                }
                break;
            }
            case 'binary':
                if (this.currentFile) {
                    return {
                        type: 'file',
                        file: this.currentFile,
                        format: this.currentFormat
                    };
                }
                break;
        }
        return null;
    }

    /**
     * Set the current file with format
     * @param {File} file - File to set
     * @param {FileFormat} format - File format
     * @returns {void}
     */
    setCurrentFile(file, format) {
        this.currentFile = file;
        this.currentFormat = format;
        this.validateFormat();
    }

    /**
     * Get the textarea ID for the current tab
     * @returns {string|null} Textarea element ID or null for binary tab
     */
    getTextareaId() {
        switch (this.activeTab) {
            case 'sam-dict':
                return 'header-input';
            case 'assembly-report':
                return 'assembly-input';
            case 'vcf':
                return 'vcf-input';
            default:
                return null;
        }
    }
}