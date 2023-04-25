function graph_link_formatter(value, row) {
  let graph_id = value;
  let link = `<a data-toggle="tooltip" data-placement="top" title="graph structure" href='../graphs/graph_${graph_id}.pdf' target="_blank" rel="noopener noreferrer">${value}</a>`;
  return link;
}
