function gene_card_link_formatter(value, row) {
    $(function () {
        $('[data-toggle="tooltip"]').tooltip()
    });
        let genes = value.split(",");
        let rows = "";
        for (let g of genes) {
            rows = `${rows}<tr><td><a data-toggle="tooltip" title="Linkout to genecards" href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=${g}" target="_blank" rel="noopener noreferrer">${g}</a></td></tr>`;
        };
        let table = `<div style="overflow-y: auto; max-height: 30vh; min-width: 140px;"><table class="table"><thead><tr><th scope="col">Genes</th></tr></thead><tbody>${rows}</tbody></table></div>`;
        let button = `<a href='#' tabindex='0' class='btn btn-primary' role='button' data-toggle='popover' data-placement="left" data-trigger='focus' data-html='true' data-content='${table}'><svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-table" viewBox="0 0 16 16">
        <path d="M0 2a2 2 0 0 1 2-2h12a2 2 0 0 1 2 2v12a2 2 0 0 1-2 2H2a2 2 0 0 1-2-2V2zm15 2h-4v3h4V4zm0 4h-4v3h4V8zm0 4h-4v3h3a1 1 0 0 0 1-1v-2zm-5 3v-3H6v3h4zm-5 0v-3H1v2a1 1 0 0 0 1 1h3zm-4-4h4V8H1v3zm0-4h4V4H1v3zm5-3v3h4V4H6zm4 4H6v3h4V8z"/>
        </svg></a>`;
        return button;
}
