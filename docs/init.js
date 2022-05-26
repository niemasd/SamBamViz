const DATA_PER_POS = 5;
const COUNT = "count";
const LOG = "log";
const PROPORTION = "proportion";
const FAS = "fasta";
const SBV = "sambamviz";

const CHECK_BOUNDS_ERR = "Please make sure bounds are valid."
const CHECK_FILES_ERR = "Please make sure file types are correct."
const CHECK_DATA_ERR = "Please make sure that file data is corect."

const BOUNDS = {
    low: 0,
    high: -1,
};

const COLORS = {
    A: '#65F73E',
    C: '#FFB340',
    G: '#EB413C',
    T: '#3C89EE',
    N: '#000000'
};

const FILE_STATUS = {
    sbv: false,
    fas: false
}

const Y_AXIS_LABEL_LOG = "Count (Log Scale)";
const Y_AXIS_LABEL_LIN = "Count (Linear SCale)";
const Y_AXIS_LABEL_PROP = "Proportion";

const TOOLTIP_COUNT = "Count";
const TOOLTIP_PROP = "Proportion";

const FAS_SELECTOR = document.getElementById('fasta-load');
const SBV_SELECTOR = document.getElementById('sambamviz-load');
const LOAD_FILES_BUTTON = document.getElementById('load-files');
const UPDATE_GRAPH_BUTTON = document.getElementById('update-graph');
const LOW_BOUNDS_SETTER = document.getElementById('lower-bound');
const UPPER_BOUNDS_SETTER = document.getElementById('upper-bound');
const GRAPH_DISPLAY = document.getElementById('view');

const COUNT_OPTION = document.getElementById('linear-option');
const LOG_OPTION = document.getElementById('log-option');
const PROP_OPTION = document.getElementById('proportion-option')

const COLOR_PICKERS = document.querySelectorAll('.color-picker');
const COLOR_PICKER_A = document.getElementById('a-color');
const COLOR_PICKER_C = document.getElementById('c-color');
const COLOR_PICKER_G = document.getElementById('g-color');
const COLOR_PICKER_T = document.getElementById('t-color');
const COLOR_PICKER_N = document.getElementById('n-color');

let tsvString;
let genSeq = [];
let tsvDataCount = [];
let tsvDataProp = [];