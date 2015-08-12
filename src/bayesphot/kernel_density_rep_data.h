/*********************************************************************
Copyright (C) 2015 Mark Krumholz
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/

/* This is just a data file that holds a tabulated function used by
   the kernel_density_rep_code. The function stored is the logarithm
   of the normalized maximum error that results from approximating a
   pair of Gaussians of equal dispersion with relative weights w and
   (1-w) and peaks separated by a distance dx (measured in units of
   the dispersion) and approximating them by a single Gaussian with
   the same mean and variance. The error is defined as the maximum
   of the absolute value of the difference between the single and
   double Gaussians over all points x, normalized to the maximum of
   the single Gaussian approximation.

   This table was generated by MRK using mathematica on 8/10/15.
*/

const double wmin = 0.0;  /* Minimum w in the table; actually used
			     10^-4, since w=0 exactly gives an error
			     of 0 exactly */
const double wmax = 0.5;  /* Maximum w in the table */
const int nw = 21;        /* Number of entries in the w direction on
			     the table */
const double dw = 0.025;  /* Size of cells in w */
const double logdxmin = -2.0; /* Minimum value of log dx in the table
				 */
const double logdxmax = 1.0; /* Maximum value of log dx in the table
				*/
const int ndx = 51;        /* Number of dx values in the table */
const double dlogdx = 0.06;  /* Size of cells in log dx */

/* The main data table */
const double gauss_err_tab[21][51] =
{{-8.30108,  -8.18109,  -8.06109,  -7.9411,  -7.82111,  -7.70112,  -7.58113,  -7.46115,  -7.34117,  -7.2212,  -7.10125,  -6.9813,  -6.86137,  -6.74147,  -6.62159,  -6.50176,  -6.38198,  -6.26226,  -6.14264,  -6.02314,  -5.9038,  -5.78466,  -5.6658,  -5.5473,  -5.42928,  -5.31188,  -5.1953,  -5.0798,  -4.9657,  -4.85343,  -4.74356,  -4.63678,  -4.53399,  -4.43629,  -4.34502,  -4.26176,  -4.05433,  -3.95661,  -3.86516,  -3.7813,  -3.70619,  -3.64039,  -3.58311,  -3.53183,  -3.48306,  -3.43417,  -3.38416,  -3.33291,  -3.28051,  -3.22707,  -3.17264},
{-5.91409,  -5.7941,  -5.6741,  -5.55411,  -5.43411,  -5.31412,  -5.19414,  -5.07415,  -4.95417,  -4.8342,  -4.71424,  -4.59428,  -4.47435,  -4.35443,  -4.23454,  -4.11469,  -3.99488,  -3.87513,  -3.75546,  -3.6359,  -3.51484,  -3.39508,  -3.27539,  -3.15581,  -3.03636,  -2.91708,  -2.79804,  -2.67929,  -2.56094,  -2.44311,  -2.32597,  -2.20972,  -2.09465,  -1.98111,  -1.86955,  -1.76057,  -1.65488,  -1.55335,  -1.45701,  -1.3669,  -1.28392,  -1.20833,  -1.13911,  -1.07346,  -1.00749,  -0.938053,  -0.863721,  -0.784172,  -0.699482,  -0.609907,  -0.515836},
{-5.62434,  -5.50435,  -5.38435,  -5.26436,  -5.14436,  -5.02437,  -4.90438,  -4.78439,  -4.66441,  -4.54444,  -4.42447,  -4.30451,  -4.18456,  -4.06464,  -3.94473,  -3.82486,  -3.70502,  -3.58524,  -3.46469,  -3.3448,  -3.22494,  -3.10514,  -2.98539,  -2.86573,  -2.74617,  -2.62675,  -2.50752,  -2.38853,  -2.26986,  -2.15161,  -2.03392,  -1.91695,  -1.80094,  -1.68618,  -1.57304,  -1.46202,  -1.35372,  -1.2489,  -1.14843,  -1.05324,  -0.964151,  -0.881443,  -0.80424,  -0.72999,  -0.65502,  -0.576337,  -0.492731,  -0.404228,  -0.31129,  -0.214516,  -0.114567},
{-5.45983,  -5.33984,  -5.21984,  -5.09984,  -4.97985,  -4.85986,  -4.73986,  -4.61988,  -4.49989,  -4.37991,  -4.25994,  -4.13998,  -4.02002,  -3.90008,  -3.78017,  -3.66027,  -3.54042,  -3.42003,  -3.3001,  -3.18018,  -3.06029,  -2.94044,  -2.82064,  -2.7009,  -2.58124,  -2.46169,  -2.34228,  -2.22306,  -2.1041,  -1.98546,  -1.86725,  -1.74962,  -1.63274,  -1.51685,  -1.40228,  -1.28942,  -1.17881,  -1.07112,  -0.967129,  -0.867743,  -0.773785,  -0.685658,  -0.602739,  -0.522822,  -0.442511,  -0.358934,  -0.270976,  -0.178826,  -0.0831315,  0.0153463,  0.115799},
{-5.34679,  -5.2268,  -5.1068,  -4.9868,  -4.86681,  -4.74681,  -4.62682,  -4.50683,  -4.38684,  -4.26686,  -4.14688,  -4.02691,  -3.90695,  -3.787,  -3.66707,  -3.54716,  -3.42728,  -3.30693,  -3.18698,  -3.06704,  -2.94712,  -2.82723,  -2.70737,  -2.58756,  -2.46781,  -2.34813,  -2.22856,  -2.10913,  -1.98988,  -1.87087,  -1.75219,  -1.63393,  -1.51625,  -1.39934,  -1.28346,  -1.16895,  -1.05628,  -0.946068,  -0.83906,  -0.736146,  -0.638209,  -0.545807,  -0.45859,  -0.374706,  -0.291029,  -0.204758,  -0.114787,  -0.0213831,  0.0746751,  0.172485,  0.271083},
{-5.26212,  -5.14212,  -5.02212,  -4.90212,  -4.78213,  -4.66213,  -4.54214,  -4.42215,  -4.30216,  -4.18217,  -4.06219,  -3.94221,  -3.82225,  -3.70229,  -3.58235,  -3.46242,  -3.34218,  -3.22221,  -3.10223,  -2.98227,  -2.86233,  -2.74239,  -2.62248,  -2.5026,  -2.38276,  -2.26297,  -2.14325,  -2.02361,  -1.9041,  -1.78475,  -1.66563,  -1.5468,  -1.42838,  -1.31052,  -1.19345,  -1.07745,  -0.962953,  -0.850521,  -0.740896,  -0.63498,  -0.533731,  -0.437865,  -0.347305,  -0.260539,  -0.174706,  -0.0870574,  0.00355283,  0.0968121,  0.191825,  0.287541,  0.382821},
{-5.19552,  -5.07552,  -4.95553,  -4.83553,  -4.71553,  -4.59554,  -4.47554,  -4.35555,  -4.23556,  -4.11557,  -3.99558,  -3.8756,  -3.75563,  -3.63566,  -3.51571,  -3.39577,  -3.27555,  -3.15556,  -3.03558,  -2.91559,  -2.79562,  -2.67565,  -2.55569,  -2.43575,  -2.31582,  -2.19592,  -2.07606,  -1.95624,  -1.83648,  -1.71682,  -1.59728,  -1.47792,  -1.35882,  -1.24011,  -1.12195,  -1.00461,  -0.888477,  -0.774088,  -0.662182,  -0.553689,  -0.449644,  -0.350919,  -0.257691,  -0.168785,  -0.081606,  0.00654293,  0.0968672,  0.189033,  0.282031,  0.374658,  0.465615},
{-5.14154,  -5.02154,  -4.90154,  -4.78154,  -4.66155,  -4.54155,  -4.42155,  -4.30156,  -4.18157,  -4.06157,  -3.94159,  -3.8216,  -3.70162,  -3.58165,  -3.46168,  -3.34154,  -3.22153,  -3.10153,  -2.98153,  -2.86153,  -2.74153,  -2.62153,  -2.50152,  -2.38152,  -2.26152,  -2.14152,  -2.02152,  -1.90153,  -1.78155,  -1.66159,  -1.54167,  -1.42183,  -1.30211,  -1.18261,  -1.06347,  -0.944924,  -0.827328,  -0.71121,  -0.597309,  -0.486588,  -0.380159,  -0.279039,  -0.183647,  -0.0931304,  -0.0051782,  0.0828407,  0.172203,  0.262574,  0.35284,  0.441665,  0.527637},
{-5.09691,  -4.97691,  -4.85691,  -4.73691,  -4.61692,  -4.49692,  -4.37692,  -4.25693,  -4.13693,  -4.01694,  -3.89694,  -3.77696,  -3.65697,  -3.53699,  -3.41701,  -3.29688,  -3.17687,  -3.05686,  -2.93685,  -2.81683,  -2.6968,  -2.57677,  -2.45673,  -2.33667,  -2.21659,  -2.0965,  -1.97638,  -1.85623,  -1.73604,  -1.61581,  -1.49555,  -1.37526,  -1.25497,  -1.13476,  -1.01473,  -0.895104,  -0.776205,  -0.658557,  -0.542908,  -0.430253,  -0.321777,  -0.218632,  -0.12146,  -0.0297286,  0.0585766,  0.145994,  0.233877,  0.321914,  0.408917,  0.493459,  0.574123},
{-5.05955,  -4.93955,  -4.81955,  -4.69955,  -4.57955,  -4.45955,  -4.33955,  -4.21956,  -4.09956,  -3.97956,  -3.85957,  -3.73957,  -3.61958,  -3.4996,  -3.37961,  -3.2595,  -3.13948,  -3.01946,  -2.89943,  -2.77939,  -2.65935,  -2.53928,  -2.4192,  -2.2991,  -2.17896,  -2.05878,  -1.93854,  -1.81825,  -1.69787,  -1.57739,  -1.45681,  -1.33612,  -1.21532,  -1.09446,  -0.973638,  -0.853037,  -0.73298,  -0.613981,  -0.496798,  -0.382461,  -0.27222,  -0.16735,  -0.0687039,  0.023942,  0.112286,  0.198743,  0.284744,  0.370037,  0.453404,  0.533413,  0.608799},
{-5.02803,  -4.90803,  -4.78803,  -4.66803,  -4.54803,  -4.42803,  -4.30803,  -4.18803,  -4.06803,  -3.94804,  -3.82804,  -3.70804,  -3.58805,  -3.46805,  -3.34806,  -3.22796,  -3.10793,  -2.98791,  -2.86787,  -2.74781,  -2.62775,  -2.50766,  -2.38754,  -2.26739,  -2.14719,  -2.02693,  -1.9066,  -1.78617,  -1.66562,  -1.54492,  -1.42405,  -1.30299,  -1.18173,  -1.06029,  -0.938752,  -0.817287,  -0.696203,  -0.576014,  -0.457489,  -0.341688,  -0.229924,  -0.123581,  -0.0237063,  0.0696222,  0.157772,  0.242999,  0.32681,  0.409061,  0.488578,  0.564041,  0.634537},
{-5.00136,  -4.88136,  -4.76136,  -4.64136,  -4.52136,  -4.40136,  -4.28136,  -4.16136,  -4.04136,  -3.92136,  -3.80136,  -3.68136,  -3.56136,  -3.44136,  -3.32136,  -3.20127,  -3.08124,  -2.9612,  -2.84115,  -2.72109,  -2.601,  -2.48089,  -2.36074,  -2.24055,  -2.1203,  -1.99997,  -1.87954,  -1.75899,  -1.63828,  -1.51738,  -1.39626,  -1.27487,  -1.15319,  -1.03124,  -0.909066,  -0.786835,  -0.664847,  -0.543615,  -0.423919,  -0.306849,  -0.193771,  -0.0861689,  0.0147307,  0.108561,  0.196348,  0.280152,  0.361553,  0.440573,  0.516185,  0.587332,  0.653638},
{-4.97881,  -4.85881,  -4.73881,  -4.61881,  -4.49881,  -4.37881,  -4.25881,  -4.13881,  -4.01881,  -3.8988,  -3.7788,  -3.6588,  -3.5388,  -3.41879,  -3.29878,  -3.1787,  -3.05867,  -2.93862,  -2.81857,  -2.69849,  -2.57839,  -2.45825,  -2.33808,  -2.21785,  -2.09755,  -1.97715,  -1.85664,  -1.73599,  -1.61514,  -1.49406,  -1.37271,  -1.25103,  -1.12898,  -1.00657,  -0.883842,  -0.76094,  -0.638162,  -0.516021,  -0.39531,  -0.277145,  -0.162939,  -0.0542642,  0.0474891,  0.141681,  0.228989,  0.311249,  0.390106,  0.465821,  0.537637,  0.604917,  0.667933},
{-4.95984,  -4.83984,  -4.71984,  -4.59984,  -4.47984,  -4.35984,  -4.23984,  -4.11984,  -3.99983,  -3.87983,  -3.75983,  -3.63982,  -3.51981,  -3.3998,  -3.27979,  -3.15972,  -3.03968,  -2.91963,  -2.79956,  -2.67948,  -2.55936,  -2.43921,  -2.31901,  -2.19874,  -2.0784,  -1.95795,  -1.83737,  -1.71662,  -1.59565,  -1.47442,  -1.35286,  -1.23092,  -1.10856,  -0.985749,  -0.862533,  -0.739048,  -0.615586,  -0.492663,  -0.371082,  -0.251979,  -0.136812,  -0.027229,  0.075232,  0.169677,  0.256441,  0.337104,  0.413373,  0.485821,  0.554109,  0.618142,  0.678806},
{-4.94405,  -4.82405,  -4.70405,  -4.58405,  -4.46404,  -4.34404,  -4.22404,  -4.10404,  -3.98403,  -3.86403,  -3.74402,  -3.62402,  -3.50401,  -3.38399,  -3.26398,  -3.14391,  -3.02387,  -2.90381,  -2.78374,  -2.66364,  -2.54351,  -2.42334,  -2.30312,  -2.18283,  -2.06245,  -1.94196,  -1.82131,  -1.70048,  -1.57941,  -1.45804,  -1.3363,  -1.21415,  -1.09151,  -0.968355,  -0.844723,  -0.72074,  -0.596697,  -0.473109,  -0.35079,  -0.230897,  -0.11492,  -0.00457799,  0.098464,  0.19308,  0.279275,  0.35836,  0.43209,  0.501431,  0.56659,  0.62808,  0.687219},
{-4.93112,  -4.81112,  -4.69112,  -4.57111,  -4.45111,  -4.33111,  -4.21111,  -4.09111,  -3.9711,  -3.8511,  -3.73109,  -3.61108,  -3.49107,  -3.37105,  -3.25103,  -3.13097,  -3.01093,  -2.89087,  -2.77078,  -2.65068,  -2.53054,  -2.41036,  -2.29012,  -2.1698,  -2.04939,  -1.92886,  -1.80816,  -1.68726,  -1.56609,  -1.44461,  -1.32273,  -1.20039,  -1.07752,  -0.954079,  -0.830098,  -0.705699,  -0.58117,  -0.45703,  -0.334099,  -0.213551,  -0.096906,  0.0140604,  0.117572,  0.212299,  0.29794,  0.375532,  0.446872,  0.513378,  0.575908,  0.635547,  0.693789},
{-4.92082,  -4.80082,  -4.68081,  -4.56081,  -4.44081,  -4.32081,  -4.20081,  -4.0808,  -3.9608,  -3.84079,  -3.72078,  -3.60077,  -3.48076,  -3.36074,  -3.24071,  -3.12068,  -3.00061,  -2.88055,  -2.76046,  -2.64035,  -2.5202,  -2.40001,  -2.27976,  -2.15942,  -2.03899,  -1.91842,  -1.79768,  -1.67672,  -1.55548,  -1.4339,  -1.31191,  -1.18941,  -1.06635,  -0.942681,  -0.818417,  -0.693682,  -0.56876,  -0.444174,  -0.320751,  -0.199676,  -0.0824955,  0.0289698,  0.132851,  0.227645,  0.312781,  0.389036,  0.458235,  0.522283,  0.582731,  0.641118,  0.698883},
{-4.91297,  -4.79297,  -4.67297,  -4.55297,  -4.43297,  -4.31296,  -4.19296,  -4.07296,  -3.95295,  -3.83294,  -3.71294,  -3.59292,  -3.47291,  -3.35289,  -3.23286,  -3.11282,  -2.99276,  -2.87269,  -2.7526,  -2.63249,  -2.51233,  -2.39213,  -2.27186,  -2.15151,  -2.03106,  -1.91046,  -1.78969,  -1.66869,  -1.54739,  -1.42574,  -1.30365,  -1.18104,  -1.05784,  -0.933984,  -0.809502,  -0.684507,  -0.559283,  -0.434354,  -0.310553,  -0.189074,  -0.0714827,  0.0403635,  0.144523,  0.239356,  0.324066,  0.399202,  0.466607,  0.528661,  0.58757,  0.645173,  0.702707},
{-4.90745,  -4.78745,  -4.66745,  -4.54745,  -4.42745,  -4.30744,  -4.18744,  -4.06744,  -3.94743,  -3.82742,  -3.70741,  -3.5874,  -3.46738,  -3.34736,  -3.22733,  -3.10729,  -2.98724,  -2.86716,  -2.74707,  -2.62695,  -2.50679,  -2.38658,  -2.26631,  -2.14595,  -2.02548,  -1.90487,  -1.78407,  -1.66303,  -1.5417,  -1.42,  -1.29785,  -1.17515,  -1.05184,  -0.927859,  -0.803222,  -0.678042,  -0.552604,  -0.427432,  -0.303363,  -0.181598,  -0.0637175,  0.0483971,  0.152752,  0.247604,  0.331992,  0.406282,  0.472333,  0.532921,  0.590791,  0.647944,  0.705379},
{-4.90417,  -4.78417,  -4.66417,  -4.54417,  -4.42417,  -4.30417,  -4.18416,  -4.06416,  -3.94415,  -3.82414,  -3.70413,  -3.58412,  -3.4641,  -3.34408,  -3.22405,  -3.10401,  -2.98395,  -2.86388,  -2.74379,  -2.62367,  -2.5035,  -2.38329,  -2.26301,  -2.14265,  -2.02217,  -1.90154,  -1.78073,  -1.65968,  -1.53832,  -1.41659,  -1.29439,  -1.17165,  -1.04827,  -0.924218,  -0.799488,  -0.674198,  -0.548632,  -0.423315,  -0.299086,  -0.177151,  -0.0590979,  0.0531762,  0.157646,  0.252507,  0.336694,  0.410458,  0.475664,  0.535358,  0.592636,  0.649563,  0.706961},
 {-4.90309,  -4.78309,  -4.66309,  -4.54308,  -4.42308,  -4.30308,  -4.18308,  -4.06307,  -3.94307,  -3.82306,  -3.70305,  -3.58303,  -3.46302,  -3.34299,  -3.22296,  -3.10292,  -2.98286,  -2.86279,  -2.7427,  -2.62258,  -2.50241,  -2.3822,  -2.26192,  -2.14155,  -2.02107,  -1.90044,  -1.77962,  -1.65856,  -1.5372,  -1.41546,  -1.29325,  -1.17049,  -1.04709,  -0.92301,  -0.798249,  -0.672922,  -0.547314,  -0.421949,  -0.297666,  -0.175675,  -0.0575645,  0.0547626,  0.15927,  0.254133,  0.338252,  0.411838,  0.476757,  0.536151,  0.593236,  0.650095,  0.707485}};
