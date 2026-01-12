/**
 * @fileoverview Split View Manager for contig comparison
 * @module managers/SplitViewManager
 */

import { escapeHtml, debounce, DEBOUNCE_DELAY } from '../utils/helpers.js';

/**
 * @typedef {'exact'|'renamed'|'conflict'|'missing'|'unknown'} MatchStatus
 * Possible match status values for contigs
 */

/**
 * @typedef {Object} ContigData
 * @property {string} name - Contig name
 * @property {number} length - Contig length in base pairs
 * @property {string} [md5] - MD5 hash of sequence
 * @property {MatchStatus} match_status - Match status
 * @property {string} [sequence_role] - Role of the sequence
 * @property {string[]} [aliases] - Alternative names
 * @property {string} [description] - Description
 * @property {Object} [assembly_info] - Assembly information
 */

/**
 * @typedef {Object} MappingEntry
 * @property {number} query_index - Index in query contigs
 * @property {number} reference_index - Index in reference contigs
 * @property {string} [query_name] - Query contig name (for renamed)
 * @property {string} [reference_name] - Reference contig name (for renamed)
 */

/**
 * @typedef {Object} MappingsData
 * @property {MappingEntry[]} exact_matches - Exact match mappings
 * @property {MappingEntry[]} renamed_matches - Renamed match mappings
 * @property {MappingEntry[]} conflicts - Conflict mappings
 * @property {number[]} query_only - Indices of query-only contigs
 * @property {number[]} reference_only - Indices of reference-only contigs
 */

/**
 * @typedef {Object} DetailedQueryData
 * @property {ContigData[]} contigs - Query contigs array
 */

/**
 * @typedef {Object} DetailedReferenceData
 * @property {string} display_name - Reference display name
 * @property {ContigData[]} contigs - Reference contigs array
 */

/**
 * @typedef {Object} DetailedComparisonData
 * @property {DetailedQueryData} query - Query data
 * @property {DetailedReferenceData} reference - Reference data
 * @property {MappingsData} mappings - Mapping information
 */

/**
 * @typedef {Object} OriginalInput
 * @property {'text'|'file'} type - Input type
 * @property {string} [content] - Text content (for text type)
 * @property {File} [file] - File object (for file type)
 * @property {string} [filename] - Filename
 */

/**
 * @typedef {'exact'|'renamed'|'conflict'} ConnectionType
 * Types of connections between contigs
 */

/**
 * Split View Manager
 * Manages the contig comparison modal with split view between query and reference contigs
 * @class
 */
export class SplitViewManager {
    /**
     * Creates a new SplitViewManager instance
     * @constructor
     */
    constructor() {
        /** @type {number|null} Index of current match being compared */
        this.currentMatchIndex = null;

        /** @type {DetailedComparisonData|null} Current detailed comparison data */
        this.currentData = null;

        /** @type {number} Current page for query contigs */
        this.queryPage = 0;

        /** @type {number} Current page for reference contigs */
        this.refPage = 0;

        /** @type {number} Number of contigs per page */
        this.pageSize = 50;

        /** @type {string} Search term for query contigs */
        this.querySearchTerm = '';

        /** @type {string} Search term for reference contigs */
        this.refSearchTerm = '';

        /** @type {ContigData[]} Filtered query contigs */
        this.filteredQueryContigs = [];

        /** @type {ContigData[]} Filtered reference contigs */
        this.filteredRefContigs = [];

        /** @type {Object|null} Original results for detailed requests */
        this.originalResults = null;

        /** @type {OriginalInput|null} Original input data */
        this.originalInput = null;

        /** @type {Object|null} Original config data */
        this.originalConfig = null;

        // Create debounced versions of filter functions for search inputs
        /** @type {function(string): void} Debounced query filter */
        this.debouncedFilterQuery = debounce(
            (term) => this.filterQueryContigs(term),
            DEBOUNCE_DELAY
        );

        /** @type {function(string): void} Debounced reference filter */
        this.debouncedFilterRef = debounce(
            (term) => this.filterReferenceContigs(term),
            DEBOUNCE_DELAY
        );
    }

    /**
     * Open the contig comparison modal for a specific match
     * @param {number} matchIndex - Index of the match to compare
     * @returns {Promise<void>}
     */
    async openComparison(matchIndex) {
        // Force refresh the page if we don't have original results
        if (!this.originalResults || !this.originalInput || !this.originalConfig) {
            // Show error in a user-friendly way
            const errorDiv = document.createElement('div');
            errorDiv.className = 'error';
            errorDiv.style.cssText = 'position: fixed; top: 20px; right: 20px; background: var(--error); color: white; padding: 1rem; border-radius: 6px; z-index: 10000;';
            errorDiv.textContent = 'Please refresh the page and run Identify again before using Compare Contigs.';
            document.body.appendChild(errorDiv);

            // Auto-remove after 5 seconds
            setTimeout(() => {
                if (errorDiv.parentNode) {
                    errorDiv.parentNode.removeChild(errorDiv);
                }
            }, 5000);
            return;
        }

        this.currentMatchIndex = matchIndex;

        // Show modal
        document.getElementById('split-view-modal').style.display = 'block';

        // Show loading state
        document.getElementById('modal-title').textContent = 'Loading contig comparison...';

        try {
            // Reconstruct the FormData the same way as the original request
            const formData = new FormData();

            if (this.originalInput.type === 'text') {
                formData.append('header_text', this.originalInput.content);
                if (this.originalInput.filename) {
                    formData.append('filename', this.originalInput.filename);
                }
            } else if (this.originalInput.type === 'file') {
                formData.append('file', this.originalInput.file);
            }

            formData.append('config', JSON.stringify(this.originalConfig));

            // Make API call to get detailed contig breakdown
            const response = await fetch(`/api/identify?mode=detailed&match_id=${matchIndex}`, {
                method: 'POST',
                body: formData
            });

            if (!response.ok) {
                const errorText = await response.text();
                throw new Error(`Failed to get detailed data: ${response.status} - ${errorText}`);
            }

            try {
                this.currentData = await response.json();
            } catch (jsonError) {
                throw new Error(`Invalid JSON response from server: ${jsonError.message}`);
            }

            this.renderSplitView();
        } catch (error) {
            document.getElementById('modal-title').textContent = 'Error loading contig comparison';
            document.getElementById('query-contigs').innerHTML = `<div class="error">Failed to load data: ${escapeHtml(error.message)}</div>`;
        }
    }

    /**
     * Render the split view comparison interface
     * @returns {void}
     * @throws {Error} If data is invalid or missing
     */
    renderSplitView() {
        const data = this.currentData;

        // Validate data structure
        if (!data) {
            throw new Error('No currentData available for split view');
        }
        if (!data.query || !data.query.contigs) {
            throw new Error('Invalid query data structure');
        }
        if (!data.reference || !data.reference.display_name) {
            throw new Error('Invalid reference data structure');
        }

        // Update modal title
        // textContent is safe - automatically escapes HTML
        document.getElementById('modal-title').textContent =
            `Contig Comparison - ${data.reference.display_name}`;

        document.getElementById('reference-title').textContent =
            `Reference Contigs (${data.reference.display_name})`;

        // Apply initial filters
        this.applyFilters();

        // Render contigs and connections
        this.renderQueryContigs();
        this.renderReferenceContigs();
        this.renderConnections();
    }

    /**
     * Apply search filters to query and reference contigs
     * @returns {void}
     */
    applyFilters() {
        const data = this.currentData;

        // Validate data before filtering
        if (!data || !data.query || !data.query.contigs) {
            this.filteredQueryContigs = [];
        } else {
            // Filter query contigs
            this.filteredQueryContigs = data.query.contigs.filter(contig =>
                contig && contig.name && contig.name.toLowerCase().includes(this.querySearchTerm.toLowerCase())
            );
        }

        if (!data || !data.reference || !data.reference.contigs) {
            this.filteredRefContigs = [];
        } else {
            // Filter reference contigs
            this.filteredRefContigs = data.reference.contigs.filter(contig =>
                contig && contig.name && contig.name.toLowerCase().includes(this.refSearchTerm.toLowerCase())
            );
        }
    }

    /**
     * Render query contigs list with pagination
     * @returns {void}
     */
    renderQueryContigs() {
        // Ensure filteredQueryContigs is an array
        if (!Array.isArray(this.filteredQueryContigs)) {
            this.filteredQueryContigs = [];
        }

        const startIdx = this.queryPage * this.pageSize;
        const endIdx = Math.min(startIdx + this.pageSize, this.filteredQueryContigs.length);
        const pageContigs = this.filteredQueryContigs.slice(startIdx, endIdx);

        let html = '';
        for (let i = 0; i < pageContigs.length; i++) {
            const contig = pageContigs[i];
            const statusClass = this.getStatusClass(contig.match_status);
            const contigIndex = startIdx + i;

            html += `
                <div class="contig-item ${statusClass}" data-contig-index="${contigIndex}" data-contig-type="query">
                    <div class="contig-name">${escapeHtml(contig.name)}</div>
                    <div class="contig-details">
                        <span class="contig-length">${contig.length.toLocaleString()} bp</span>
                        <span class="contig-status ${statusClass}">${escapeHtml(contig.match_status)}</span>
                        ${contig.md5 ? `<span class="contig-md5" title="${escapeHtml(contig.md5)}">${contig.md5.substring(0, 8)}...</span>` : ''}
                        <span class="click-hint">Click for details</span>
                    </div>
                </div>
            `;
        }

        document.getElementById('query-contigs').innerHTML = html;

        // Add event listeners for query contig clicks
        document.querySelectorAll('#query-contigs .contig-item').forEach((item) => {
            item.addEventListener('click', () => {
                const contigIndex = parseInt(item.dataset.contigIndex, 10);
                const contig = this.filteredQueryContigs[contigIndex];
                this.showContigDetails('query', contig);
            });
        });

        this.updateQueryPagination();
    }

    /**
     * Render reference contigs list with pagination
     * @returns {void}
     */
    renderReferenceContigs() {
        const startIdx = this.refPage * this.pageSize;
        const endIdx = Math.min(startIdx + this.pageSize, this.filteredRefContigs.length);
        const pageContigs = this.filteredRefContigs.slice(startIdx, endIdx);

        let html = '';
        for (let i = 0; i < pageContigs.length; i++) {
            const contig = pageContigs[i];
            const statusClass = this.getStatusClass(contig.match_status);
            const contigIndex = startIdx + i;

            html += `
                <div class="contig-item ${statusClass}" data-contig-index="${contigIndex}" data-contig-type="reference">
                    <div class="contig-name">${escapeHtml(contig.name)}</div>
                    <div class="contig-details">
                        <span class="contig-length">${contig.length.toLocaleString()} bp</span>
                        <span class="contig-status ${statusClass}">${escapeHtml(contig.match_status)}</span>
                        ${contig.md5 ? `<span class="contig-md5" title="${escapeHtml(contig.md5)}">${contig.md5.substring(0, 8)}...</span>` : ''}
                        <span class="click-hint">Click for details</span>
                    </div>
                </div>
            `;
        }

        document.getElementById('reference-contigs').innerHTML = html;

        // Add event listeners for reference contig clicks
        document.querySelectorAll('#reference-contigs .contig-item').forEach((item) => {
            item.addEventListener('click', () => {
                const contigIndex = parseInt(item.dataset.contigIndex, 10);
                const contig = this.filteredRefContigs[contigIndex];
                this.showContigDetails('reference', contig);
            });
        });

        this.updateReferencePagination();
    }

    /**
     * Render SVG connection lines between matched contigs
     * @returns {void}
     */
    renderConnections() {
        const svg = document.getElementById('connections-svg');

        // Clear existing connections
        svg.innerHTML = '';

        const data = this.currentData;
        const queryStartIdx = this.queryPage * this.pageSize;
        const refStartIdx = this.refPage * this.pageSize;

        // Draw connections for exact matches
        data.mappings.exact_matches.forEach(match => {
            if (this.isInCurrentPage(match.query_index, queryStartIdx, this.pageSize) &&
                this.isInCurrentPage(match.reference_index, refStartIdx, this.pageSize)) {
                this.drawConnection(svg, match.query_index - queryStartIdx, match.reference_index - refStartIdx, 'exact');
            }
        });

        // Draw connections for renamed matches
        data.mappings.renamed_matches.forEach(match => {
            if (this.isInCurrentPage(match.query_index, queryStartIdx, this.pageSize) &&
                this.isInCurrentPage(match.reference_index, refStartIdx, this.pageSize)) {
                this.drawConnection(svg, match.query_index - queryStartIdx, match.reference_index - refStartIdx, 'renamed');
            }
        });

        // Draw connections for conflicts
        data.mappings.conflicts.forEach(match => {
            if (this.isInCurrentPage(match.query_index, queryStartIdx, this.pageSize) &&
                match.reference_index && this.isInCurrentPage(match.reference_index, refStartIdx, this.pageSize)) {
                this.drawConnection(svg, match.query_index - queryStartIdx, match.reference_index - refStartIdx, 'conflict');
            }
        });
    }

    /**
     * Check if an index is within the current page range
     * @param {number} index - Index to check
     * @param {number} startIdx - Start index of current page
     * @param {number} pageSize - Page size
     * @returns {boolean} True if index is in current page
     */
    isInCurrentPage(index, startIdx, pageSize) {
        return index >= startIdx && index < startIdx + pageSize;
    }

    /**
     * Draw a connection line between two contigs
     * @param {SVGSVGElement} svg - SVG element to draw in
     * @param {number} queryIdx - Query contig index (relative to page)
     * @param {number} refIdx - Reference contig index (relative to page)
     * @param {ConnectionType} type - Type of connection
     * @returns {void}
     */
    drawConnection(svg, queryIdx, refIdx, type) {
        const line = document.createElementNS('http://www.w3.org/2000/svg', 'path');

        // Calculate positions (simplified - would need actual element positions)
        const queryY = queryIdx * 60 + 30; // Approximate contig item height
        const refY = refIdx * 60 + 30;
        const svgWidth = svg.clientWidth;
        const startX = 0;
        const endX = svgWidth;
        const midX = svgWidth / 2;

        // Create curved path
        const path = `M ${startX} ${queryY} C ${midX} ${queryY} ${midX} ${refY} ${endX} ${refY}`;

        line.setAttribute('d', path);
        line.setAttribute('stroke', this.getConnectionColor(type));
        line.setAttribute('stroke-width', '2');
        line.setAttribute('fill', 'none');
        line.setAttribute('class', `connection connection-${type}`);

        svg.appendChild(line);
    }

    /**
     * Get CSS color variable for connection type
     * @param {ConnectionType} type - Connection type
     * @returns {string} CSS variable string
     */
    getConnectionColor(type) {
        switch (type) {
            case 'exact': return 'var(--success)';
            case 'renamed': return 'var(--accent)';
            case 'conflict': return 'var(--error)';
            default: return 'var(--text-muted)';
        }
    }

    /**
     * Get CSS class for match status
     * @param {MatchStatus} status - Match status
     * @returns {string} CSS class name
     */
    getStatusClass(status) {
        switch (status) {
            case 'exact': return 'status-exact';
            case 'renamed': return 'status-renamed';
            case 'conflict': return 'status-conflict';
            case 'missing': return 'status-missing';
            default: return 'status-unknown';
        }
    }

    /**
     * Update pagination info and controls for query contigs
     * @returns {void}
     */
    updateQueryPagination() {
        const totalPages = Math.ceil(this.filteredQueryContigs.length / this.pageSize);
        const currentPage = this.queryPage + 1;

        document.getElementById('query-pagination-info').textContent =
            `Page ${currentPage} of ${totalPages} (${this.filteredQueryContigs.length} contigs)`;

        // Add pagination controls
        this.renderPaginationControls('query-pagination', this.queryPage, totalPages,
            (page) => {
                this.queryPage = page;
                this.renderQueryContigs();
                this.renderConnections();
            });
    }

    /**
     * Update pagination info and controls for reference contigs
     * @returns {void}
     */
    updateReferencePagination() {
        const totalPages = Math.ceil(this.filteredRefContigs.length / this.pageSize);
        const currentPage = this.refPage + 1;

        document.getElementById('ref-pagination-info').textContent =
            `Page ${currentPage} of ${totalPages} (${this.filteredRefContigs.length} contigs)`;

        // Add pagination controls
        this.renderPaginationControls('ref-pagination', this.refPage, totalPages,
            (page) => {
                this.refPage = page;
                this.renderReferenceContigs();
                this.renderConnections();
            });
    }

    /**
     * Render pagination control buttons
     * @param {string} containerId - ID of container element
     * @param {number} currentPage - Current page number (0-indexed)
     * @param {number} totalPages - Total number of pages
     * @param {function(number): void} onPageChange - Callback when page changes
     * @returns {void}
     */
    renderPaginationControls(containerId, currentPage, totalPages, onPageChange) {
        const container = document.getElementById(containerId);

        if (totalPages <= 1) {
            container.innerHTML = '';
            return;
        }

        // Clear container
        container.innerHTML = '';

        // Previous button
        if (currentPage > 0) {
            const prevBtn = document.createElement('button');
            prevBtn.className = 'page-btn';
            prevBtn.textContent = 'Prev';
            prevBtn.onclick = () => onPageChange(currentPage - 1);
            container.appendChild(prevBtn);
        }

        // Page numbers (show max 5 around current)
        const start = Math.max(0, currentPage - 2);
        const end = Math.min(totalPages, currentPage + 3);

        for (let i = start; i < end; i++) {
            const pageBtn = document.createElement('button');
            pageBtn.className = i === currentPage ? 'page-btn active' : 'page-btn';
            pageBtn.textContent = i + 1;
            pageBtn.onclick = () => onPageChange(i);
            container.appendChild(pageBtn);
        }

        // Next button
        if (currentPage < totalPages - 1) {
            const nextBtn = document.createElement('button');
            nextBtn.className = 'page-btn';
            nextBtn.textContent = 'Next';
            nextBtn.onclick = () => onPageChange(currentPage + 1);
            container.appendChild(nextBtn);
        }
    }

    /**
     * Filter query contigs by search term
     * @param {string} searchTerm - Search term to filter by
     * @returns {void}
     */
    filterQueryContigs(searchTerm) {
        this.querySearchTerm = searchTerm;
        this.queryPage = 0; // Reset to first page
        this.applyFilters();
        this.renderQueryContigs();
        this.renderConnections();
    }

    /**
     * Filter reference contigs by search term
     * @param {string} searchTerm - Search term to filter by
     * @returns {void}
     */
    filterReferenceContigs(searchTerm) {
        this.refSearchTerm = searchTerm;
        this.refPage = 0; // Reset to first page
        this.applyFilters();
        this.renderReferenceContigs();
        this.renderConnections();
    }

    /**
     * Close the comparison modal and reset state
     * @returns {void}
     */
    closeComparison() {
        document.getElementById('split-view-modal').style.display = 'none';
        this.currentMatchIndex = null;
        this.currentData = null;
        this.queryPage = 0;
        this.refPage = 0;
        this.querySearchTerm = '';
        this.refSearchTerm = '';

        // Clear search inputs
        document.getElementById('query-search').value = '';
        document.getElementById('ref-search').value = '';
    }

    /**
     * Show detailed information modal for a contig
     * @param {'query'|'reference'} type - Type of contig (query or reference)
     * @param {ContigData|string} contigData - Contig data object or JSON string
     * @returns {void}
     */
    showContigDetails(type, contigData) {
        try {
            // Parse the contig data
            const contig = typeof contigData === 'string' ? JSON.parse(contigData.replace(/&quot;/g, '"')) : contigData;

            // Create details modal content
            const modalHtml = `
                <div id="contig-details-modal" class="modal" style="display: block; z-index: 1100;">
                    <div class="modal-content" style="max-width: 600px;">
                        <div class="modal-header">
                            <h2>${escapeHtml(contig.name)} Details</h2>
                            <span class="modal-close" onclick="document.getElementById('contig-details-modal').remove()">&times;</span>
                        </div>
                        <div class="modal-body" style="padding: 1.5rem;">
                            <div class="contig-detail-grid">
                                <!-- Always show name and length -->
                                <div class="detail-row">
                                    <span class="detail-label">Name:</span>
                                    <span class="detail-value">${escapeHtml(contig.name)}</span>
                                </div>
                                <div class="detail-row">
                                    <span class="detail-label">Length:</span>
                                    <span class="detail-value">${contig.length.toLocaleString()} bp</span>
                                </div>

                                <!-- Only show MD5 if it exists -->
                                ${contig.md5 ? `
                                <div class="detail-row">
                                    <span class="detail-label">MD5 Hash:</span>
                                    <span class="detail-value monospace">${escapeHtml(contig.md5)}</span>
                                </div>
                                ` : ''}

                                <!-- Only show sequence role if it exists -->
                                ${contig.sequence_role ? `
                                <div class="detail-row">
                                    <span class="detail-label">Sequence Role:</span>
                                    <span class="detail-value">${escapeHtml(contig.sequence_role)}</span>
                                </div>
                                ` : ''}

                                <!-- Only show aliases if they exist -->
                                ${contig.aliases && contig.aliases.length > 0 ? `
                                <div class="detail-row">
                                    <span class="detail-label">Aliases:</span>
                                    <span class="detail-value">${contig.aliases.map(alias => escapeHtml(alias)).join(', ')}</span>
                                </div>
                                ` : ''}

                                <!-- Only show description if it exists -->
                                ${contig.description ? `
                                <div class="detail-row">
                                    <span class="detail-label">Description:</span>
                                    <span class="detail-value">${escapeHtml(contig.description)}</span>
                                </div>
                                ` : ''}

                                <!-- Only show assembly info if it exists -->
                                ${contig.assembly_info ? `
                                <div class="detail-row">
                                    <span class="detail-label">Assembly:</span>
                                    <span class="detail-value">${escapeHtml(contig.assembly_info.assembly)}${contig.assembly_info.organism ? ` (${escapeHtml(contig.assembly_info.organism)})` : ''}</span>
                                </div>
                                ` : ''}

                                <!-- Always show match status and source -->
                                <div class="detail-row">
                                    <span class="detail-label">Match Status:</span>
                                    <span class="detail-value status-badge ${this.getStatusClass(contig.match_status)}">${escapeHtml(contig.match_status)}</span>
                                </div>
                                <div class="detail-row">
                                    <span class="detail-label">Source:</span>
                                    <span class="detail-value">${type === 'query' ? 'Your Input' : 'Reference Genome'}</span>
                                </div>
                            </div>
                            <div style="margin-top: 1.5rem; text-align: center;">
                                <button onclick="document.getElementById('contig-details-modal').remove()" style="
                                    background: var(--success);
                                    color: white;
                                    border: none;
                                    padding: 0.75rem 1.5rem;
                                    border-radius: 6px;
                                    cursor: pointer;
                                    font-size: 1rem;
                                ">Close</button>
                            </div>
                        </div>
                    </div>
                </div>
            `;

            // Add modal to page
            document.body.insertAdjacentHTML('beforeend', modalHtml);

        } catch (error) {
            // Show error in a more user-friendly way
            const errorDiv = document.createElement('div');
            errorDiv.className = 'error';
            errorDiv.style.cssText = 'position: fixed; top: 20px; right: 20px; background: var(--error); color: white; padding: 1rem; border-radius: 6px; z-index: 10000;';
            errorDiv.textContent = 'Error showing contig details: ' + error.message;
            document.body.appendChild(errorDiv);

            // Auto-remove after 5 seconds
            setTimeout(() => {
                if (errorDiv.parentNode) {
                    errorDiv.parentNode.removeChild(errorDiv);
                }
            }, 5000);
        }
    }
}