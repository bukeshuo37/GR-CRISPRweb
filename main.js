function updateEditors(){
    // read mode from radio group
    let mode = document.querySelector('input[name="mode"]:checked').value;
    let editorSelect = document.getElementById("editor");
    editorSelect.innerHTML="";
    let cas9 = ["WT-SpCas9","eSpCas9(1.1)","SpCas9-HF1"];
    let be = ["ABE7.10","ABEmax","ABE8e","BE4","CBE4max","Target-AID"];
    let list = (mode=="cas9") ? cas9 : be;
    list.forEach(e=>{
        let option=document.createElement("option");
        option.value=e; option.text=e;
        editorSelect.appendChild(option);
    });
    // Set default selected option
    if(editorSelect.options.length > 0) {
        editorSelect.selectedIndex = 0;
    }
    updateSortOptions();
}

function updateSortOptions(){
    let mode = document.querySelector('input[name="mode"]:checked').value;
    let sortSelect = document.getElementById("sort_by");
    if(!sortSelect) return;
    sortSelect.innerHTML = "";
    let opts = [];
    if(mode=="cas9"){
        opts = ["","GR-CRISPRScore","DeepMEnsScore","TransCrisprScore"];
    } else {
        opts = ["","GR-CRISPRScore","DeepBEScore","BEDeeponScore"];
    }
    opts.forEach(v=>{
        let o=document.createElement('option');
        o.value=v; o.text=v||'Sort (Default: Name → GR-CRISPRScore)';
        sortSelect.appendChild(o);
    });
}
function preloadData() {

    fetch("/api/search",{
        method:"POST",
        headers:{"Content-Type":"application/json"},
        body:JSON.stringify({omim_id:"", mutation:"", mode:"cas9", editor:"WT-SpCas9"})
    })
    .then(r=>r.json())
    .then(data=>{
        console.log("Data preloaded successfully");
    })
    .catch(error=>{
        console.error("Error preloading data:", error);
    });
}
document.addEventListener("DOMContentLoaded",function(){
    const params = new URLSearchParams(window.location.search);
    const modeParam = params.get("mode");
    if(modeParam === "cas9" || modeParam === "be") {
        const modeEl = document.querySelector(`input[name="mode"][value="${modeParam}"]`);
        if(modeEl) modeEl.checked = true;
    }

    updateEditors();
    document.querySelectorAll('input[name="mode"]').forEach(r=>r.addEventListener('change',updateEditors));

    const omimEl = document.getElementById("omim");
    const mutationEl = document.getElementById("mutation");
    const editorEl = document.getElementById("editor");
    const sortEl = document.getElementById("sort_by");

    const omimParam = params.get("omim_id");
    const mutationParam = params.get("mutation");
    const editorParam = params.get("editor");
    const sortParam = params.get("sort_by");
    const autoSearch = params.get("autosearch") === "1";

    if(omimEl && omimParam) omimEl.value = omimParam;
    if(mutationEl && mutationParam) mutationEl.value = mutationParam;
    if(editorEl && editorParam && Array.from(editorEl.options).some(o => o.value === editorParam)) {
        editorEl.value = editorParam;
    }
    if(sortEl && sortParam && Array.from(sortEl.options).some(o => o.value === sortParam)) {
        sortEl.value = sortParam;
    }

    // Initialize editors first
    updateEditors();
    
    // Check if input fields have default values after a short delay to ensure everything is rendered
    setTimeout(() => {
        const omimValue = omimEl ? omimEl.value : "";
        const mutationValue = mutationEl ? mutationEl.value : "";
        const editorValue = document.getElementById("editor").value;
        
        console.log('Page loaded - omimValue:', omimValue);
        console.log('Page loaded - mutationValue:', mutationValue);
        console.log('Page loaded - editorValue:', editorValue);
        console.log('Page loaded - autoSearch:', autoSearch);
        console.log('Page loaded - omimParam:', omimParam);
        console.log('Page loaded - mutationParam:', mutationParam);
        
        if(autoSearch || omimParam || mutationParam || omimValue || mutationValue) {
            console.log('Triggering search...');
            search();
        } else {
            console.log('Preloading data...');
            preloadData();
        }
    }, 300);
});

// Pagination variables
let currentData = [];
let allSearchData = [];
let currentPage = 1;
const itemsPerPage = 20;
let searchMode = '';
let beTargetFilter = null;

function shouldMergeColumn(columnName) {
    if(searchMode !== "be") {
        return true;
    }
    const beMergeColumns = new Set(["Name", "sgRNA", "PAM"]);
    return beMergeColumns.has(columnName);
}

function normalizeMergeValue(value) {
    if(value === undefined || value === null) {
        return "";
    }
    if(typeof value === "string") {
        return value.trim();
    }
    return String(value).trim();
}

function normalizeTargetSequence(value) {
    if(value === undefined || value === null) {
        return "";
    }
    return String(value).trim().toUpperCase().substring(0, 20);
}

function formatDisplayValue(key, value) {
    if(value === undefined || value === null) {
        return "";
    }
    if(key && key.includes("Score")) {
        const num = Number(value);
        if(Number.isFinite(num)) {
            return num.toFixed(3);
        }
    }
    return value;
}

function filterBeByTarget(nameEncoded, sgRNAEncoded) {
    if(searchMode !== "be") {
        return;
    }

    const targetName = decodeURIComponent(nameEncoded || "").trim();
    const targetSgRNA = decodeURIComponent(sgRNAEncoded || "").trim();
    

    let omim = document.getElementById("omim").value;
    let mutation = document.getElementById("mutation").value;
    let mode = document.querySelector('input[name="mode"]:checked').value;
    let editor = document.getElementById("editor").value;
    

    fetch("/api/search_original", {
        method: "POST",
        headers: {"Content-Type": "application/json"},
        body: JSON.stringify({omim_id: omim, mutation: mutation, mode: mode, editor: editor, name: targetName, sgRNA: targetSgRNA})
    })
    .then(r => {
        if (!r.ok) {
            throw new Error('Network response was not ok');
        }
        return r.json();
    })
    .then(data => {
        if (data.error) {
            alert('Error: ' + data.error);
            return;
        }
        currentData = data;
        allSearchData = data; 
        currentPage = 1;
        beTargetFilter = { name: targetName, sgRNA: targetSgRNA };
        renderPage(1);
    })
    .catch(error => {
        console.error('Error:', error);
        alert('An error occurred while fetching data. Please try again.');
    });
}

function clearBeTargetFilter() {

    search();
    beTargetFilter = null;
}

function renderPage(pageNum) {
    const totalPages = Math.ceil(currentData.length / itemsPerPage);
    if(pageNum < 1) pageNum = 1;
    if(pageNum > totalPages) pageNum = totalPages;
    currentPage = pageNum;

    const startIdx = (pageNum - 1) * itemsPerPage;
    const endIdx = Math.min(startIdx + itemsPerPage, currentData.length);
    const pageData = currentData.slice(startIdx, endIdx);

    let html = "";
    if(pageData.length === 0) {
        html = "<p style='color:#888;'>No data to display</p>";
    } else {
        if(searchMode === "be" && beTargetFilter) {
            html += "<div style='margin-bottom:10px;padding:10px;background:#eef6ff;border:1px solid #c9defc;border-radius:6px;'>";
            html += "Filtered by Name: <strong>" + beTargetFilter.name + "</strong>, sgRNA: <strong>" + beTargetFilter.sgRNA + "</strong> ";
            html += "<button onclick='clearBeTargetFilter()' style='margin-left:10px;'>Clear Filter</button>";
            html += "</div>";
        }


        let columns;
        if(beTargetFilter) {

            columns = Object.keys(pageData[0]);
        } else {

            columns = Object.keys(pageData[0]).filter(col => col !== "Outcome");
        }
        let displayColumns = columns.map(col => {
            if(searchMode === "be" && col === "sgRNA") {
                return "sgRNA";
            }
            if(searchMode === "be" && col === "Outcome") {
                return "Outcome";
            }
            return col;
        });

        // Calculate rowspan for cell merging
        let rowspanMap = {}; // rowspanMap[rowIdx][colIdx] = {span: count, skip: boolean}
        for(let i = 0; i < pageData.length; i++) {
            rowspanMap[i] = {};
            for(let j = 0; j < columns.length; j++) {
                rowspanMap[i][j] = {span: 1, skip: false};
            }
        }

        // Detect cells that should be merged
        for(let colIdx = 0; colIdx < columns.length; colIdx++) {
            if(!shouldMergeColumn(columns[colIdx])) {
                continue;
            }
            let rowIdx = 0;
            while(rowIdx < pageData.length) {
                let currentVal = normalizeMergeValue(pageData[rowIdx][columns[colIdx]]);
                let spanCount = 1;

                // Count consecutive rows with same value
                let nextRowIdx = rowIdx + 1;
                while(nextRowIdx < pageData.length &&
                        normalizeMergeValue(pageData[nextRowIdx][columns[colIdx]]) === currentVal) {
                    spanCount++;
                    nextRowIdx++;
                }

                // Set rowspan for first row and skip flag for subsequent rows
                if(spanCount > 1) {
                    rowspanMap[rowIdx][colIdx].span = spanCount;
                    for(let k = 1; k < spanCount; k++) {
                        rowspanMap[rowIdx + k][colIdx].skip = true;
                    }
                }

                rowIdx = nextRowIdx;
            }
        }

        html += "<table><tr>";
        displayColumns.forEach(col => html += "<th>" + col + "</th>");
        html += "</tr>";

        pageData.forEach((row, rowIdx) => {
            html += "<tr>";
            columns.forEach((k, colIdx) => {
                // Skip if this cell should be merged with the one above
                if(rowspanMap[rowIdx][colIdx].skip) {
                    return;
                }

                let val = formatDisplayValue(k, row[k]);

                if(k === "sgRNA") {
                    if(typeof val === "string" && val.length > 20) {
                        val = val.substring(0, 20);
                    }

                    if(searchMode === "be") {
                        const encodedName = encodeURIComponent(String(row.Name ?? "").trim());
                        const encodedSgRNA = encodeURIComponent(String(val ?? "").trim());
                        val = `<a href="#" onclick="filterBeByTarget('${encodedName}','${encodedSgRNA}');return false;" style="color: blue; text-decoration: underline;">${val}</a>`;
                    }
                }

                let rowspan = rowspanMap[rowIdx][colIdx].span;
                let rowspanAttr = rowspan > 1 ? ` rowspan="${rowspan}"` : '';
                let mergedStyle = (rowspan > 1 && shouldMergeColumn(k))
                    ? " style='text-align:center;vertical-align:middle;'"
                    : "";
                html += "<td" + rowspanAttr + mergedStyle + ">" + val + "</td>";
            });
            html += "</tr>";
        });
        html += "</table>";
    }

    // Add pagination controls AFTER table
    html += "<div class='pagination-controls' style='margin-top: 20px;'>";
    html += "<button onclick='goToPage(1)' " + (pageNum===1 ? "disabled" : "") + ">⬅ First</button>";
    html += "<button onclick='goToPage(" + (pageNum-1) + ")' " + (pageNum===1 ? "disabled" : "") + ">‹ Previous</button>";

    // Page numbers
    html += "<div class='pagination-button-group' style='margin: 0 10px;'>";
    const maxButtons = 7;
    let startPage = Math.max(1, pageNum - Math.floor(maxButtons/2));
    let endPage = Math.min(totalPages, startPage + maxButtons - 1);
    if(endPage - startPage < maxButtons - 1) {
        startPage = Math.max(1, endPage - maxButtons + 1);
    }

    if(startPage > 1) html += "<span class='pagination-ellipsis'>...</span>";
    for(let i = startPage; i <= endPage; i++) {
        if(i === pageNum) {
            html += `<span class='page-number active' onclick='goToPage(${i})'>${i}</span>`;
        } else {
            html += `<span class='page-number' onclick='goToPage(${i})'>${i}</span>`;
        }
    }
    if(endPage < totalPages) html += "<span class='pagination-ellipsis'>...</span>";
    html += "</div>";

    html += "<button onclick='goToPage(" + (pageNum+1) + ")' " + (pageNum===totalPages ? "disabled" : "") + ">Next ›</button>";
    html += "<button onclick='goToPage(" + totalPages + ")' " + (pageNum===totalPages ? "disabled" : "") + ">Last ⮕</button>";

    html += "<div style='margin-top:10px; color:#666; font-size:14px;'>";
    html += "Total: " + currentData.length + " records | Page " + pageNum + " of " + totalPages;
    html += "</div>";
    html += "</div>";

    document.getElementById("result").innerHTML = html;
}

function goToPage(pageNum) {
    renderPage(pageNum);
}

function search(){
    let omim = document.getElementById("omim").value;
    let mutation = document.getElementById("mutation").value;
    let mode = document.querySelector('input[name="mode"]:checked').value;
    let editorSelect = document.getElementById("editor");
    let editor = editorSelect.value || (editorSelect.options.length > 0 ? editorSelect.options[0].value : "WT-SpCas9");
    let sortBy = document.getElementById("sort_by") ? document.getElementById("sort_by").value : "";
    
    console.log('Search triggered with:', { omim, mutation, mode, editor });

    fetch("/api/search",{
        method:"POST",
        headers:{"Content-Type":"application/json"},
        body:JSON.stringify({omim_id:omim, mutation:mutation, mode:mode, editor:editor, sort_by:sortBy})
    })
    .then(r=>r.json())
    .then(data=>{
        allSearchData = data;
        currentData = data;
        searchMode = mode;
        currentPage = 1;
        beTargetFilter = null;

        if(data.length==0){
            document.getElementById("result").innerHTML="<p style='color:#888;'>No data found for the current search criteria</p><p style='color:#888;'>(Please wait a few minutes when loading data for the first time)</p>";
        } else {
            renderPage(1);
        }
    });
}
