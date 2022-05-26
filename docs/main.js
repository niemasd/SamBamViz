/**
 * Event listener code for clicking on buttons
 */
LOAD_FILES_BUTTON.addEventListener('click', () => {
    generateData();

    let lowerDataBound = BOUNDS.low;
    let upperDataBound = BOUNDS.high;
    loadSBV(lowerDataBound, upperDataBound, tsvDataCount, tsvDataProp);
});

SBV_SELECTOR.addEventListener('change', () => {
    processFiles(SBV);
});

FAS_SELECTOR.addEventListener('change', () => {
    processFiles(FAS);
});

UPDATE_GRAPH_BUTTON.addEventListener('click', () => {
    let lowerDataBound = BOUNDS.low;
    let upperDataBound = BOUNDS.high;

    lowerDataBound = parseInt(LOW_BOUNDS_SETTER.value);
    upperDataBound = parseInt(UPPER_BOUNDS_SETTER.value) + 1;
    loadSBV(lowerDataBound, upperDataBound, tsvDataCount, tsvDataProp);
});

/**
 * Handles color changing inputs
 */
COLOR_PICKERS.forEach((input) => {
    input.addEventListener('change', () => {
        setVegaScheme();
    });
});

/**
 * Turns the data from a TSV file to an array in memory
 */
async function generateData() {
    try {
        await tsvToArr(tsvString);
    } catch (error) {
        console.log(error);
    }
}

/**
 * Processes a file differently depending on the file type
 * @param {String} fileType - the file extension
 */
async function processFiles(fileType) {
    try {
        if (fileType === SBV) {
            const sbvFile =
                SBV_SELECTOR.files[0];

            if (sbvFile.name.split(".")[1] === "tsv") {
                FILE_STATUS.sbv = true;
                await processSamBamViz(sbvFile);
            } else {
                FILE_STATUS.sbv = false;
            }
        }

        if (fileType === FAS) {
            const fasFile = FAS_SELECTOR.files[0];

            if (fasFile.name.split(".")[1] === "fas") {
                FILE_STATUS.fas = true;
                await processFasta(fasFile);
            } else {
                FILE_STATUS.fas = false;
            }
        }
    } catch (error) {
        console.log(error);
    }
}

/**
 * Takes a SamBamViz file and puts it into a string
 * @param {File} file - the SamBamViz file as a TSV
 * @returns a Promise after the file has been read
 */
function processSamBamViz(file) {
    return new Promise((resolve) => {
        let reader = new FileReader();

        reader.onload = () => {
            tsvString = reader.result;
            resolve('done');
        };

        reader.readAsText(file);
    });
}

/**
 * Takes a fasta file and parses the genome sequence into a string
 * @param {File} file - the FASTA file
 * @returns a Promise after the file has been read
 */
function processFasta(file) {
    return new Promise((resolve) => {
        let reader = new FileReader();

        reader.onload = () => {
            genSeq = [];
            let genSeqStr = reader.result;
            for (
                let i = genSeqStr.indexOf('\n') + 1;
                i < genSeqStr.length;
                i++
            ) {
                if (genSeqStr[i] != '\n') {
                    genSeq.push(genSeqStr[i]);
                }
            }
            resolve('done');
        };

        reader.readAsText(file);
    });
}

/**
 * Gets the colors from the color selector
 */
function setVegaScheme() {
    COLORS.A = COLOR_PICKER_A.value;
    COLORS.C = COLOR_PICKER_C.value;
    COLORS.G = COLOR_PICKER_G.value;
    COLORS.T = COLOR_PICKER_T.value;
    COLORS.N = COLOR_PICKER_N.value;
}

/**
 * Checks to make sure that there are no invalid inputs before loading
 * the graph
 * @param {Number} lowBounds - the lower bounds of data inputted
 * @param {Number} highBounds - the higher bounds of data inputted
 * @returns an error message depending on the error found
 */
function checkIssues(lowBounds, highBounds) {
    const fasFile = FAS_SELECTOR.files[0];

    if (!FILE_STATUS.sbv || (fasFile && !FILE_STATUS.fas)) {
        return CHECK_FILES_ERR;
    }

    if (lowBounds > highBounds) {
        return CHECK_BOUNDS_ERR;
    }

    return "";
}

/**
 * Prepares to take the data and convert it into a graph
 * @param {Number} lowerDataBound - the lower bounds of data inputted
 * @param {Number} upperDataBound - the higher bounds of data inputted
 * @param {String[]} tsvData - count data in an array
 * @param {String[]} tsvDataProp - proportional data in an array
 * @returns Nothing
 */
function loadSBV(lowerDataBound, upperDataBound, tsvData, tsvDataProp) {
    let errorMessage = checkIssues(lowerDataBound, upperDataBound);
    if (errorMessage) {
        alert(errorMessage);
        return;
    }

    vega.scheme('strands', [
        COLORS.A,
        COLORS.C,
        COLORS.G,
        COLORS.T,
        COLORS.N
    ]);

    let countData = tsvData.slice(
        lowerDataBound * DATA_PER_POS,
        upperDataBound * DATA_PER_POS
    );

    let propData = tsvDataProp.slice(
        lowerDataBound * DATA_PER_POS,
        upperDataBound * DATA_PER_POS
    );

    if (COUNT_OPTION.checked) {
        loadBarGraph(countData, COUNT);
    } else if (LOG_OPTION.checked) {
        loadBarGraph(countData, LOG);
    } else if (PROP_OPTION.checked) {
        loadBarGraph(propData, PROPORTION);
    }
}