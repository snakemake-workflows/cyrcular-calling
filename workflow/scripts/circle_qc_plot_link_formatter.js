function circle_qc_plot_link_formatter(value, row) {
  let graph_id = row["graph_id"]
  let circle_id = row["circle_id"]
  let link = `<a data-toggle="tooltip" data-placement="top" title="circle candidate quality control plot: binned coverage across circle, total and broken down by type of base alignment " href='../qc_plots/graph_${graph_id}_${circle_id}.html' target="_blank" rel="noopener noreferrer">${value}</a>`;
  return link;
}
