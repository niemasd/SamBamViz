/**
 * Loads the stacked bar graph from data extracted from TSV file and user 
 * customizations (e.g. type of graph) by using Vega. data passed into function
 * should already be in user-defined data range. Colors should already be set.
 * @param {Object[]} data - data extracted from TSV file
 * @param {String} yAxisData - type of graph
 */
function loadBarGraph(data, yAxisData) {
    let viewElement = GRAPH_DISPLAY;
    let width = 2 * viewElement.offsetWidth / 3;
    let height = viewElement.offsetHeight / 2;
    let yscale;                     // to tell vega type of graph

    let yAxisName = 'default';      // for y-axis label
    if (yAxisData == LOG) {
        yscale = 'symlog';
        yAxisName = Y_AXIS_LABEL_LOG;
    } else {
        yscale = 'linear';
        yAxisName = Y_AXIS_LABEL_LIN;
    }

    let label = 'default';          // for tooltip 
    if (yAxisData == PROPORTION) {
        label = TOOLTIP_PROP;
        yAxisName = Y_AXIS_LABEL_PROP;
    } else {
        label = TOOLTIP_COUNT;
    }

    let hasGenSeq = genSeq.length > 0;  // if reference file is available
    let tooltipDisplay = hasGenSeq ?
        `{'Letter': datum.type, 'Position': datum.pos, ${label}: datum.value, 'Reference': datum.ref}` :
        `{'Letter': datum.type, 'Position': datum.pos, ${label}: datum.value}`;

    // stacked bar chart
    var barChart = {
        $schema: 'https://vega.github.io/schema/vega/v5.json',
        description: 'A stacked bar chart for SamBamViz.',
        width: width,
        height: height,
        padding: 5,

        data: [
            {
                name: 'table',
                values: [...data],
                transform: [
                    {
                        type: 'stack',
                        groupby: ['pos'],
                        sort: { field: 'type' },
                        field: 'value',
                        sort: { field: 'value', order: "ascending" }
                    }
                ]
            }
        ],

        scales: [
            {
                name: 'x',
                type: 'band',
                range: 'width',
                domain: { data: 'table', field: 'pos' }
            },
            {
                name: 'y',
                type: yscale,
                range: 'height',
                nice: true,
                zero: true,
                domain: { data: 'table', field: 'value' }
            },
            {
                name: 'color',
                type: 'ordinal',
                range: { scheme: "strands" },
                domain: { data: 'table', field: 'type' },
            }
        ],

        legends: [
            {
                fill: 'color',
                orient: 'right',
                title: 'Codes',
                format: '',
                encode: {
                    symbols: {
                        update: {
                            shape: { value: 'square' },
                            stroke: { value: '#ccc' },
                            strokeWidth: { value: 0.2 }
                        }
                    }
                }
            }
        ],

        axes: [
            { orient: 'bottom', scale: 'x', zindex: 1, labelOverlap: "parity", title: "Genome Position" },
            { orient: 'left', scale: 'y', zindex: 1, labelOverlap: "greedy", title: `${yAxisName}` }
        ],

        marks: [
            {
                type: 'rect',
                from: { data: 'table' },
                encode: {
                    enter: {
                        x: { scale: 'x', field: 'pos' },
                        width: { scale: 'x', band: 1, offset: -1 },
                        y: { scale: 'y', field: 'y0' },
                        y2: { scale: 'y', field: 'y1' },
                        fill: { scale: 'color', field: 'type' },
                        tooltip: { "signal": tooltipDisplay }
                    },
                    update: {
                        fillOpacity: { value: 1 }
                    },
                    hover: {
                        fillOpacity: { value: 0.5 }
                    }
                }
            }
        ]
    };

    vegaEmbed('#view', barChart);
}