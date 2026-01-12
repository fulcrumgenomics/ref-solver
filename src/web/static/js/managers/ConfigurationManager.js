/**
 * @fileoverview Configuration Manager for application settings
 * @module managers/ConfigurationManager
 */

import { clamp } from '../utils/helpers.js';

/**
 * @typedef {Object} ScoringWeights
 * @property {number} md5Jaccard - Weight for MD5 Jaccard score (0-100)
 * @property {number} nameLength - Weight for name/length Jaccard score (0-100)
 * @property {number} md5Coverage - Weight for MD5 coverage score (0-100)
 * @property {number} orderScore - Weight for order score (0-100)
 */

/**
 * @typedef {Object} AppConfig
 * @property {number} scoreThreshold - Minimum score threshold (0-1)
 * @property {number} resultLimit - Maximum number of results to display
 * @property {ScoringWeights} scoringWeights - Scoring weight configuration
 */

/**
 * Configuration Manager
 * Handles application configuration storage, validation, and UI updates
 * @class
 */
export class ConfigurationManager {
    /**
     * Creates a new ConfigurationManager instance
     * @constructor
     */
    constructor() {
        /** @type {number} Score threshold (0-100) */
        this.scoreThreshold = 80;

        /** @type {number} Maximum results to display */
        this.resultLimit = 10;

        /** @type {ScoringWeights} Scoring weights */
        this.scoringWeights = {
            md5Jaccard: 40,
            nameLength: 30,
            md5Coverage: 20,
            orderScore: 10
        };

        this.loadConfig();
    }

    /**
     * Load configuration from localStorage with validation
     * @returns {void}
     */
    loadConfig() {
        // Load and validate score threshold (0-100)
        const storedThreshold = localStorage.getItem('scoreThreshold');
        if (storedThreshold !== null) {
            const parsed = parseInt(storedThreshold, 10);
            this.scoreThreshold = isNaN(parsed) ? 80 : clamp(parsed, 0, 100);
        }

        // Load and validate result limit (1-50)
        const storedLimit = localStorage.getItem('resultLimit');
        if (storedLimit !== null) {
            const parsed = parseInt(storedLimit, 10);
            this.resultLimit = isNaN(parsed) ? 10 : clamp(parsed, 1, 50);
        }

        // Load and validate scoring weights (0-100 each)
        this.scoringWeights = {
            md5Jaccard: this._loadWeight('md5JaccardWeight', 40),
            nameLength: this._loadWeight('nameLengthWeight', 30),
            md5Coverage: this._loadWeight('md5CoverageWeight', 20),
            orderScore: this._loadWeight('orderScoreWeight', 10)
        };

        this.updateUI();
    }

    /**
     * Load and validate a single weight value from localStorage
     * @private
     * @param {string} key - localStorage key
     * @param {number} defaultValue - Default value if not found or invalid
     * @returns {number} Validated weight value (0-100)
     */
    _loadWeight(key, defaultValue) {
        const stored = localStorage.getItem(key);
        if (stored === null) return defaultValue;
        const parsed = parseInt(stored, 10);
        return isNaN(parsed) ? defaultValue : clamp(parsed, 0, 100);
    }

    /**
     * Save current configuration to localStorage
     * @returns {void}
     */
    saveConfig() {
        localStorage.setItem('scoreThreshold', String(this.scoreThreshold));
        localStorage.setItem('resultLimit', String(this.resultLimit));
        localStorage.setItem('md5JaccardWeight', String(this.scoringWeights.md5Jaccard));
        localStorage.setItem('nameLengthWeight', String(this.scoringWeights.nameLength));
        localStorage.setItem('md5CoverageWeight', String(this.scoringWeights.md5Coverage));
        localStorage.setItem('orderScoreWeight', String(this.scoringWeights.orderScore));
    }

    /**
     * Update score threshold with validation
     * @param {number|string} value - New threshold value
     * @returns {void}
     */
    setScoreThreshold(value) {
        const parsed = typeof value === 'string' ? parseInt(value, 10) : value;
        this.scoreThreshold = isNaN(parsed) ? this.scoreThreshold : clamp(parsed, 0, 100);
        this.saveConfig();
    }

    /**
     * Update a scoring weight with validation
     * @param {keyof ScoringWeights} key - Weight key to update
     * @param {number|string} value - New weight value
     * @returns {void}
     */
    setWeight(key, value) {
        const parsed = typeof value === 'string' ? parseInt(value, 10) : value;
        if (!isNaN(parsed) && key in this.scoringWeights) {
            this.scoringWeights[key] = clamp(parsed, 0, 100);
            this.saveConfig();
        }
    }

    /**
     * Update UI elements to reflect current configuration
     * @returns {void}
     */
    updateUI() {
        const elements = {
            'score-threshold': this.scoreThreshold,
            'threshold-value': this.scoreThreshold + '%',
            'result-limit': this.resultLimit,
            'md5-jaccard-weight': this.scoringWeights.md5Jaccard,
            'md5-jaccard-value': this.scoringWeights.md5Jaccard + '%',
            'name-length-weight': this.scoringWeights.nameLength,
            'name-length-value': this.scoringWeights.nameLength + '%',
            'md5-coverage-weight': this.scoringWeights.md5Coverage,
            'md5-coverage-value': this.scoringWeights.md5Coverage + '%',
            'order-score-weight': this.scoringWeights.orderScore,
            'order-score-value': this.scoringWeights.orderScore + '%'
        };

        for (const [id, value] of Object.entries(elements)) {
            const el = document.getElementById(id);
            if (el) {
                if (el.tagName === 'INPUT' || el.tagName === 'SELECT') {
                    el.value = value;
                } else {
                    el.textContent = value;
                }
            }
        }
    }

    /**
     * Get configuration object for API requests
     * @returns {AppConfig} Configuration object with normalized values
     */
    getConfig() {
        return {
            scoreThreshold: this.scoreThreshold / 100,
            resultLimit: this.resultLimit,
            scoringWeights: { ...this.scoringWeights }
        };
    }

    /**
     * Reset configuration to defaults
     * @returns {void}
     */
    resetToDefaults() {
        this.scoreThreshold = 80;
        this.resultLimit = 10;
        this.scoringWeights = {
            md5Jaccard: 40,
            nameLength: 30,
            md5Coverage: 20,
            orderScore: 10
        };
        this.saveConfig();
        this.updateUI();
    }
}