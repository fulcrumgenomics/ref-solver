/**
 * @fileoverview Main Application Entry Point
 * Coordinates all managers and sets up global event handlers
 * @module main
 */

import { ConfigurationManager } from './managers/ConfigurationManager.js';
import { TabManager } from './managers/TabManager.js';
import { ResultsManager, HelpSystem } from './managers/ResultsManager.js';
import { SplitViewManager } from './managers/SplitViewManager.js';
import {
    escapeHtml,
    examples,
    validateFileSize,
    clamp,
    MAX_FILE_SIZE,
    MAX_TEXT_FILE_SIZE
} from './utils/helpers.js';

/**
 * @typedef {'text'|'assembly'|'vcf'|'binary'} FileFormatType
 * Supported file format types for upload
 */

/**
 * @typedef {'md5-jaccard'|'name-length'|'md5-coverage'|'order-score'} WeightType
 * Weight type identifiers used in UI
 */

// Create global manager instances
/** @type {ConfigurationManager} */
const configManager = new ConfigurationManager();

/** @type {TabManager} */
const tabManager = new TabManager();

/** @type {ResultsManager} */
const resultsManager = new ResultsManager();

/** @type {HelpSystem} */
const helpSystem = new HelpSystem();

/** @type {SplitViewManager} */
const splitViewManager = new SplitViewManager();

// Expose managers to window for inline event handlers
window.configManager = configManager;
window.tabManager = tabManager;
window.resultsManager = resultsManager;
window.helpSystem = helpSystem;
window.splitViewManager = splitViewManager;

// Event handler functions (exposed globally for inline handlers)

/**
 * Switch to a different input format tab
 * @param {string} tabId - ID of the tab to switch to
 * @returns {void}
 */
function showTab(tabId) {
    tabManager.switchTab(tabId);
}

/**
 * Toggle the help panel visibility
 * @returns {void}
 */
function toggleHelp() {
    helpSystem.toggle();
}

/**
 * Show a specific help tab
 * @param {string} tabId - ID of the help tab to show
 * @returns {void}
 */
function showHelpTab(tabId) {
    helpSystem.showTab(tabId);
}

/**
 * Update the score threshold value
 * @param {string|number} value - New threshold value (0-100)
 * @returns {void}
 */
function updateThreshold(value) {
    const parsed = parseInt(value, 10);
    const validated = isNaN(parsed) ? configManager.scoreThreshold : clamp(parsed, 0, 100);
    configManager.scoreThreshold = validated;
    document.getElementById('threshold-value').textContent = validated + '%';
    configManager.saveConfig();
}

/**
 * Update a scoring weight value
 * @param {WeightType} type - Weight type identifier
 * @param {string|number} value - New weight value (0-100)
 * @returns {void}
 */
function updateWeight(type, value) {
    const mapping = {
        'md5-jaccard': 'md5Jaccard',
        'name-length': 'nameLength',
        'md5-coverage': 'md5Coverage',
        'order-score': 'orderScore'
    };

    const key = mapping[type];
    if (!key) return;

    const parsed = parseInt(value, 10);
    const validated = isNaN(parsed) ? configManager.scoringWeights[key] : clamp(parsed, 0, 100);
    configManager.scoringWeights[key] = validated;
    document.getElementById(type + '-value').textContent = validated + '%';
    configManager.saveConfig();
}

/**
 * Toggle the advanced options panel visibility
 * @returns {void}
 */
function toggleAdvanced() {
    const options = document.getElementById('advanced-options');
    const arrow = document.getElementById('advanced-arrow');

    if (options.classList.contains('expanded')) {
        options.classList.remove('expanded');
        arrow.textContent = '▶';
    } else {
        options.classList.add('expanded');
        arrow.textContent = '▼';
    }
}

/**
 * Toggle expansion of a result row details
 * @param {number} index - Index of the result to toggle
 * @returns {void}
 */
function toggleResultDetails(index) {
    resultsManager.toggleExpansion(index);
}

/**
 * Toggle the input section collapse state
 * @returns {void}
 */
function toggleInputSection() {
    const inputSection = document.getElementById('input-section');
    const arrow = document.getElementById('input-collapse-arrow');

    if (inputSection.classList.contains('collapsed')) {
        inputSection.classList.remove('collapsed');
        arrow.textContent = '▼';
    } else {
        inputSection.classList.add('collapsed');
        arrow.textContent = '▲';
    }
}

/**
 * Collapse the input section after running identify
 * @returns {void}
 */
function collapseInputAfterIdentify() {
    const inputSection = document.getElementById('input-section');
    const arrow = document.getElementById('input-collapse-arrow');

    inputSection.classList.add('collapsed');
    arrow.textContent = '▲';
}

/**
 * Show an error message to the user
 * @param {string} message - Error message to display
 * @returns {void}
 */
function showUploadError(message) {
    const errorDiv = document.createElement('div');
    errorDiv.className = 'error';
    errorDiv.style.cssText = 'position: fixed; top: 20px; right: 20px; background: var(--error); color: white; padding: 1rem; border-radius: 6px; z-index: 10000; max-width: 400px;';
    errorDiv.textContent = message;
    document.body.appendChild(errorDiv);

    setTimeout(() => {
        if (errorDiv.parentNode) {
            errorDiv.parentNode.removeChild(errorDiv);
        }
    }, 5000);
}

/**
 * Handle file upload and preview
 * @param {HTMLInputElement} input - File input element
 * @param {FileFormatType} format - Expected file format
 * @returns {void}
 */
function handleFileUpload(input, format) {
    const file = input.files[0];
    if (!file) return;

    // Validate file size based on format
    const maxSize = format === 'binary' ? MAX_FILE_SIZE : MAX_TEXT_FILE_SIZE;
    const validation = validateFileSize(file, maxSize);

    if (!validation.valid) {
        showUploadError(validation.error);
        input.value = ''; // Clear the file input
        return;
    }

    tabManager.currentFile = file;
    tabManager.currentFormat = format;

    // Show preview for text-based formats
    if (['text', 'assembly', 'vcf'].includes(format)) {
        const reader = new FileReader();
        reader.onload = function(e) {
            const content = e.target.result;

            if (format === 'text') {
                document.getElementById('header-input').value = content;
                tabManager.validateFormat();
            } else if (format === 'assembly') {
                document.getElementById('assembly-input').value = content;
                tabManager.validateFormat();
            } else if (format === 'vcf') {
                document.getElementById('vcf-input').value = content;
                tabManager.validateFormat();
            }
        };
        reader.readAsText(file);
    } else if (format === 'binary') {
        document.getElementById('binary-preview').style.display = 'block';
        document.getElementById('binary-filename').textContent = file.name;
        tabManager.validateFormat();
    }
}

/**
 * Load example data into the current tab's textarea
 * @param {string} tabType - Tab type identifier (sam-dict, assembly-report, vcf)
 * @returns {void}
 */
function loadExample(tabType) {
    const example = examples[tabType];
    if (!example) {
        return;
    }

    let textareaId;
    switch (tabType) {
        case 'sam-dict':
            textareaId = 'header-input';
            break;
        case 'assembly-report':
            textareaId = 'assembly-input';
            break;
        case 'vcf':
            textareaId = 'vcf-input';
            break;
        default:
            return;
    }

    const textarea = document.getElementById(textareaId);
    if (textarea) {
        textarea.value = example;
        tabManager.validateFormat();
    }
}

/**
 * Clear a textarea and revalidate the form
 * @param {string} textareaId - ID of the textarea to clear
 * @returns {void}
 */
function clearTextArea(textareaId) {
    const textarea = document.getElementById(textareaId);
    if (textarea) {
        textarea.value = '';
        tabManager.validateFormat();
    }
}

// Expose all functions globally for inline event handlers
window.showTab = showTab;
window.toggleHelp = toggleHelp;
window.showHelpTab = showHelpTab;
window.updateThreshold = updateThreshold;
window.updateWeight = updateWeight;
window.toggleAdvanced = toggleAdvanced;
window.toggleResultDetails = toggleResultDetails;
window.toggleInputSection = toggleInputSection;
window.collapseInputAfterIdentify = collapseInputAfterIdentify;
window.handleFileUpload = handleFileUpload;
window.loadExample = loadExample;
window.clearTextArea = clearTextArea;
window.escapeHtml = escapeHtml;

// Initialize on page load
document.addEventListener('DOMContentLoaded', function() {
    configManager.updateUI();
    tabManager.validateFormat();

    // Add input event listeners for real-time validation
    const textareas = ['header-input', 'assembly-input', 'vcf-input'];
    textareas.forEach(id => {
        const textarea = document.getElementById(id);
        if (textarea) {
            textarea.addEventListener('input', function() {
                tabManager.validateFormat();
            });
            textarea.addEventListener('paste', function() {
                // Validate after paste content is processed
                setTimeout(() => tabManager.validateFormat(), 10);
            });
        }
    });

    // Add change listener for result limit
    const resultLimitEl = document.getElementById('result-limit');
    if (resultLimitEl) {
        resultLimitEl.addEventListener('change', function() {
            configManager.resultLimit = parseInt(this.value, 10);
            configManager.saveConfig();
        });
    }

    // Form submission
    const identifyForm = document.getElementById('identify-form');
    if (identifyForm) {
        identifyForm.addEventListener('submit', async function(e) {
            e.preventDefault();

            // Collapse input section after clicking identify
            collapseInputAfterIdentify();

            const input = tabManager.getCurrentInput();

            if (!input) {
                const resultsDiv = document.getElementById('results');
                resultsDiv.innerHTML = '<div class="error">No valid input found. Please enter text or upload a file.</div>';
                return;
            }

            const resultsDiv = document.getElementById('results');
            resultsDiv.innerHTML = '<div class="loading">Analyzing...</div>';

            try {
                const formData = new FormData();

                if (input.type === 'text') {
                    formData.append('header_text', input.content);
                    if (input.filename) {
                        formData.append('filename', input.filename);
                    }
                } else if (input.type === 'file') {
                    formData.append('file', input.file);
                }

                // Add configuration
                const config = configManager.getConfig();

                formData.append('config', JSON.stringify(config));
                const response = await fetch('/api/identify', {
                    method: 'POST',
                    body: formData
                });

                if (!response.ok) {
                    const responseText = await response.text();
                    throw new Error(`HTTP ${response.status}: ${response.statusText} - ${responseText}`);
                }

                let data;
                try {
                    data = await response.json();
                } catch (jsonError) {
                    throw new Error(`Invalid JSON response from server: ${jsonError.message}`);
                }

                // Store input and config in SplitViewManager for detailed requests
                splitViewManager.originalInput = input;
                splitViewManager.originalConfig = config;

                resultsManager.renderResults(data, config);
            } catch (error) {
                resultsDiv.innerHTML = `<div class="error">Error: ${escapeHtml(error.message)}</div>`;
            }
        });
    }
});

// Export for testing
export {
    configManager,
    tabManager,
    resultsManager,
    helpSystem,
    splitViewManager,
    showTab,
    toggleHelp,
    showHelpTab,
    updateThreshold,
    updateWeight,
    toggleAdvanced,
    toggleResultDetails,
    toggleInputSection,
    collapseInputAfterIdentify,
    handleFileUpload,
    loadExample,
    clearTextArea
};