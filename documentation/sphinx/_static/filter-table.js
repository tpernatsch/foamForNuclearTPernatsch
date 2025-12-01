function filterTable(input, table) {
    const value = input.value.toLowerCase();

    for (let row of table.rows) {
        if (row.rowIndex === 0) continue;  // skip header

        // Regular text inside the row
        let text = row.innerText.toLowerCase();

        // Add all <img alt="..."> values
        let alts = Array.from(row.querySelectorAll("img"))
                        .map(img => img.alt ? img.alt.toLowerCase() : "")
                        .join(" ");

        // Combined search content
        let searchable = text + " " + alts;

        row.style.display = searchable.includes(value) ? "" : "none";
    }
}

window.addEventListener("DOMContentLoaded", () => {
    const filters = document.querySelectorAll("input.tableFilter");
    const tables  = document.querySelectorAll("table.filterable-table");

    filters.forEach((input, i) => {
        const table = tables[i];  // pair input with table by order
        input.addEventListener("input", () => filterTable(input, table));
    });
});
