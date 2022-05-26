/**
 * This function takes in a string with data read from the tsv file that the 
 * user uploaded and creates two Object[], tsvDataCount and tsvDataProp, which
 * Vega will use to plot the bar charts. This function is called after pressing
 * the "Load SamBamViz" button. Clicking this button will generate the default
 * ranges on the web page.
 * @param {String} tsvData - read from uploaded tsv file as string
 */
async function tsvToArr(tsvData) {
    tsvDataCount = [];
    tsvDataProp = [];
    BOUNDS.low = 0;
    BOUNDS.high = -1;

    await d3.tsvParse(tsvData, function (data) {
        const total = parseInt(data.A) + parseInt(data.G) + parseInt(data.C) +
            parseInt(data.T) + parseInt(data.Other) || 1;
        const ref = genSeq[BOUNDS.high + 1]

        // Count Data
        tsvDataCount.push({
            value: parseInt(data.A),
            type: "A",
            pos: parseInt(data.Pos),
            ref: ref
        });
        tsvDataCount.push({
            value: parseInt(data.C),
            type: "C",
            pos: parseInt(data.Pos),
            ref: ref
        });
        tsvDataCount.push({
            value: parseInt(data.G),
            type: "G",
            pos: parseInt(data.Pos),
            ref: ref
        });
        tsvDataCount.push({
            value: parseInt(data.T),
            type: "T",
            pos: parseInt(data.Pos),
            ref: ref
        });
        tsvDataCount.push({
            value: parseInt(data.Other),
            type: "Other",
            pos: parseInt(data.Pos),
            ref: ref
        });

        // Proportion Data
        tsvDataProp.push({
            value: parseInt(data.A) / total,
            type: "A",
            pos: parseInt(data.Pos),
            ref: ref
        });
        tsvDataProp.push({
            value: parseInt(data.C) / total,
            type: "C",
            pos: parseInt(data.Pos),
            ref: ref
        });
        tsvDataProp.push({
            value: parseInt(data.G) / total,
            type: "G",
            pos: parseInt(data.Pos),
            ref: ref
        });
        tsvDataProp.push({
            value: parseInt(data.T) / total,
            type: "T",
            pos: parseInt(data.Pos),
            ref: ref
        });
        tsvDataProp.push({
            value: parseInt(data.Other) / total,
            type: "Other",
            pos: parseInt(data.Pos),
            ref: ref
        });

        BOUNDS.high++;
    });

    // update ranges on page
    LOW_BOUNDS_SETTER.value = BOUNDS.low;
    UPPER_BOUNDS_SETTER.value = BOUNDS.high;
    LOW_BOUNDS_SETTER.max = BOUNDS.high;
    UPPER_BOUNDS_SETTER.max = BOUNDS.high;
}