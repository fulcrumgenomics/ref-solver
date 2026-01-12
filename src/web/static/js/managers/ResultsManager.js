/**
 * @fileoverview Results Manager for rendering and interacting with identification results
 * @module managers/ResultsManager
 */

import { escapeHtml } from '../utils/helpers.js';

/**
 * @typedef {Object} ScoreInfo
 * @property {number} composite - Composite score (0-1)
 * @property {string} confidence - Confidence level (High, Medium, Low)
 */

/**
 * @typedef {Object} ReferenceInfo
 * @property {string} display_name - Display name of the reference
 * @property {string} assembly - Assembly name
 * @property {string} source - Source organization
 * @property {string} [download_url] - Optional download URL
 */

/**
 * @typedef {Object} Suggestion
 * @property {'rename'|'reorder'|'replace'|'use_as_is'|'realign'} type - Suggestion type
 * @property {string} [command] - Command to execute
 * @property {string} [contig] - Contig name
 * @property {string} [reason] - Reason for suggestion
 */

/**
 * @typedef {Object} Match
 * @property {ScoreInfo} score - Score information
 * @property {ReferenceInfo} reference - Reference information
 * @property {string} match_type - Type of match
 * @property {number} exact_matches - Number of exact matches
 * @property {number} renamed_matches - Number of renamed matches
 * @property {number} conflicts - Number of conflicts
 * @property {boolean} reordered - Whether contigs are reordered
 * @property {Suggestion[]} suggestions - Array of suggestions
 */

/**
 * @typedef {Object} QueryInfo
 * @property {number} contig_count - Number of contigs
 * @property {number} md5_coverage - MD5 coverage percentage (0-1)
 * @property {string} naming_convention - Detected naming convention
 */

/**
 * @typedef {Object} ResultsData
 * @property {QueryInfo} query - Query information
 * @property {Match[]} matches - Array of matches
 * @property {string} [error] - Error message if any
 */

/**
 * @typedef {Object} AppConfig
 * @property {number} scoreThreshold - Score threshold (0-1)
 * @property {number} resultLimit - Maximum results to display
 */

/**
 * Results Manager
 * Handles rendering and interaction with identification results
 * @class
 */
export class ResultsManager {
    /**
     * Creates a new ResultsManager instance
     * @constructor
     */
    constructor() {
        /** @type {Match[]} Current results array */
        this.currentResults = [];

        /** @type {Set<number>} Set of expanded item indices */
        this.expandedItems = new Set();
    }

    /**
     * Render identification results to the DOM
     * @param {ResultsData} data - Results data from API
     * @param {AppConfig} config - Application configuration
     * @returns {void}
     */
    renderResults(data, config) {
        const container = document.getElementById('results');
        this.currentResults = data.matches || [];

        // Store results in SplitViewManager for detailed requests
        window.splitViewManager.originalResults = data;

        if (data.error) {
            container.innerHTML = `<div class="error">${escapeHtml(data.error)}</div>`;
            return;
        }

        // Filter by score threshold
        const filteredMatches = this.currentResults.filter(match =>
            match.score.composite >= config.scoreThreshold
        );

        let html = `
            <div class="stats">
                <div>Contigs: <span>${data.query.contig_count}</span></div>
                <div>MD5 Coverage: <span>${(data.query.md5_coverage * 100).toFixed(0)}%</span></div>
                <div>Naming: <span>${escapeHtml(data.query.naming_convention)}</span></div>
                <div>Matches: <span>${filteredMatches.length}</span></div>
            </div>
        `;

        if (filteredMatches.length === 0) {
            html += '<div class="error">No matches found above the score threshold</div>';
        } else {
            html += this.renderExpandableTable(filteredMatches);
        }

        container.innerHTML = html;
    }

    /**
     * Render the expandable results table
     * @param {Match[]} matches - Array of filtered matches
     * @returns {string} HTML string for the table
     */
    renderExpandableTable(matches) {
        let html = `
            <div class="results-summary">
                <div class="expand-controls">
                    <h3 style="color: var(--accent); margin: 0;">Results Table</h3>
                    <button class="expand-all-btn" onclick="resultsManager.toggleExpandAll()">
                        <span id="expand-all-text">Expand All</span>
                    </button>
                </div>
                <table class="summary-table">
                    <thead>
                        <tr>
                            <th>Reference</th>
                            <th>Assembly</th>
                            <th>Source</th>
                            <th>Score</th>
                            <th>Confidence</th>
                            <th>Match Type</th>
                        </tr>
                    </thead>
                    <tbody>
        `;

        for (let i = 0; i < matches.length; i++) {
            const match = matches[i];
            const confidence = match.score.confidence.toLowerCase();
            const isExpanded = this.expandedItems.has(i);

            // Main row
            html += `
                <tr onclick="resultsManager.toggleExpansion(${i})" data-result-index="${i}" class="${isExpanded ? 'selected' : ''}">
                    <td><strong>${escapeHtml(match.reference.display_name)}</strong></td>
                    <td>${escapeHtml(match.reference.assembly)}</td>
                    <td>${escapeHtml(match.reference.source)}</td>
                    <td>${(match.score.composite * 100).toFixed(1)}%</td>
                    <td><span class="confidence-badge confidence-${escapeHtml(confidence)}">${escapeHtml(match.score.confidence)}</span></td>
                    <td>${escapeHtml(match.match_type)}</td>
                </tr>
            `;

            // Expandable row with details
            html += `
                <tr class="expandable-row ${isExpanded ? 'expanded' : ''}" id="expandable-${i}">
                    <td colspan="6">
                        <div class="expandable-content">
                            <div class="result-meta">
                                <div class="meta-item">Assembly: <span>${escapeHtml(match.reference.assembly)}</span></div>
                                <div class="meta-item">Source: <span>${escapeHtml(match.reference.source)}</span></div>
                                <div class="meta-item">Match Type: <span>${escapeHtml(match.match_type)}</span></div>
                                <div class="meta-item">Score: <span>${(match.score.composite * 100).toFixed(1)}%</span></div>
                                <div class="meta-item">Exact Matches: <span>${match.exact_matches}</span></div>
                                <div class="meta-item">Renamed: <span>${match.renamed_matches}</span></div>
                                ${match.conflicts > 0 ? `<div class="meta-item" style="color: var(--error)">Conflicts: <span>${match.conflicts}</span></div>` : ''}
                                ${match.reordered ? `<div class="meta-item" style="color: var(--warning)">Reordered: <span>Yes</span></div>` : ''}
                            </div>
                            ${this.renderSuggestions(match.suggestions)}
                            <div style="margin-top: 1rem;">
                                <button class="compare-button" onclick="splitViewManager.openComparison(${i})" style="
                                    background: var(--success);
                                    color: white;
                                    border: none;
                                    padding: 0.5rem 1rem;
                                    border-radius: 6px;
                                    cursor: pointer;
                                    font-size: 0.9rem;
                                    transition: background 0.2s ease;
                                " onmouseover="this.style.background='#27ae60'" onmouseout="this.style.background='var(--success)'">
                                    Compare Contigs
                                </button>
                            </div>
                            ${match.reference.download_url ? `<div style="margin-top: 1rem; font-size: 0.85rem;"><a href="${escapeHtml(match.reference.download_url)}" target="_blank" style="color: var(--accent);">Download Reference</a></div>` : ''}
                        </div>
                    </td>
                </tr>
            `;
        }

        html += '</tbody></table></div>';
        return html;
    }

    /**
     * Render suggestions section for a match
     * @param {Suggestion[]} suggestions - Array of suggestions
     * @returns {string} HTML string for suggestions
     */
    renderSuggestions(suggestions) {
        if (!suggestions || suggestions.length === 0) return '';

        let html = '<div class="suggestions"><h4>Suggestions</h4>';

        for (const s of suggestions) {
            switch (s.type) {
                case 'rename':
                    html += `<p><strong>Rename contigs:</strong></p><div class="suggestion-code">${escapeHtml(s.command)}</div>`;
                    break;
                case 'reorder':
                    html += `<p><strong>Reorder contigs:</strong></p><div class="suggestion-code">${escapeHtml(s.command)}</div>`;
                    break;
                case 'replace':
                    html += `<p><strong>Replace ${escapeHtml(s.contig)}:</strong> ${escapeHtml(s.reason)}</p>`;
                    break;
                case 'use_as_is':
                    html += `<p style="color: var(--success)"><strong>Safe to use as-is</strong></p>`;
                    break;
                case 'realign':
                    html += `<p style="color: var(--error)"><strong>Realignment needed:</strong> ${escapeHtml(s.reason)}</p>`;
                    break;
            }
        }

        html += '</div>';
        return html;
    }

    /**
     * Toggle expansion state of a result row
     * @param {number} index - Index of the result to toggle
     * @returns {void}
     */
    toggleExpansion(index) {
        if (this.expandedItems.has(index)) {
            this.expandedItems.delete(index);
        } else {
            this.expandedItems.add(index);
        }

        // Update UI
        const expandableRow = document.getElementById(`expandable-${index}`);
        const summaryRow = document.querySelector(`tr[data-result-index="${index}"]`);

        if (this.expandedItems.has(index)) {
            expandableRow?.classList.add('expanded');
            summaryRow?.classList.add('selected');
        } else {
            expandableRow?.classList.remove('expanded');
            summaryRow?.classList.remove('selected');
        }

        this.updateExpandAllButton();
    }

    /**
     * Toggle expansion state of all result rows
     * @returns {void}
     */
    toggleExpandAll() {
        const allExpanded = this.expandedItems.size === this.currentResults.length;

        if (allExpanded) {
            // Collapse all
            this.expandedItems.clear();
            document.querySelectorAll('.expandable-row').forEach(row => {
                row.classList.remove('expanded');
            });
            document.querySelectorAll('tr[data-result-index]').forEach(row => {
                row.classList.remove('selected');
            });
        } else {
            // Expand all
            for (let i = 0; i < this.currentResults.length; i++) {
                this.expandedItems.add(i);
            }
            document.querySelectorAll('.expandable-row').forEach(row => {
                row.classList.add('expanded');
            });
            document.querySelectorAll('tr[data-result-index]').forEach(row => {
                row.classList.add('selected');
            });
        }

        this.updateExpandAllButton();
    }

    /**
     * Update the expand/collapse all button text
     * @returns {void}
     */
    updateExpandAllButton() {
        const expandAllText = document.getElementById('expand-all-text');
        if (!expandAllText) return;

        const allExpanded = this.expandedItems.size === this.currentResults.length;
        expandAllText.textContent = allExpanded ? 'Collapse All' : 'Expand All';
    }
}

/**
 * @typedef {'overview'|'formats'|'scoring'|'troubleshooting'} HelpTabId
 * Valid help tab identifiers
 */

/**
 * Help System
 * Manages the expandable help panel
 * @class
 */
export class HelpSystem {
    /**
     * Creates a new HelpSystem instance
     * @constructor
     */
    constructor() {
        /** @type {boolean} Whether help panel is expanded */
        this.isExpanded = false;

        /** @type {HelpTabId} Currently active help tab */
        this.activeTab = 'overview';
    }

    /**
     * Toggle the help panel visibility
     * @returns {void}
     */
    toggle() {
        this.isExpanded = !this.isExpanded;
        const helpContent = document.getElementById('help-content');
        const helpToggle = document.querySelector('.help-toggle');
        const arrow = document.getElementById('help-arrow');

        if (this.isExpanded) {
            helpContent.classList.add('expanded');
            helpToggle.classList.add('expanded');
            arrow.textContent = '▲';
        } else {
            helpContent.classList.remove('expanded');
            helpToggle.classList.remove('expanded');
            arrow.textContent = '▼';
        }
    }

    /**
     * Show a specific help tab
     * @param {HelpTabId} tabId - ID of the tab to show
     * @returns {void}
     */
    showTab(tabId) {
        // Remove active class from all tabs and contents
        document.querySelectorAll('.help-tab').forEach(tab => {
            tab.classList.remove('active');
        });
        document.querySelectorAll('.help-tab-content').forEach(content => {
            content.classList.remove('active');
        });

        // Add active class to selected tab and content
        event.target.classList.add('active');
        document.getElementById('help-' + tabId).classList.add('active');

        this.activeTab = tabId;
    }
}