using Toybox.WatchUi as Ui;
using Toybox.Graphics as Gfx;
using Toybox.Application as App;
using Toybox.Application.Properties as Props;
using Toybox.Application.Storage as Storage;
using Toybox.Time as Time;
using Toybox.Position as Position;
using Toybox.Time.Gregorian;
using Toybox.Math;
using Toybox.System as Sys;
using Toybox.Lang;

using Astro as astro;

(:glance)
module Phase {

//    JDN stands for Julian Day Number
//    Angles here are in degrees
//    1980 January 0.0 in JDN
    const EPOCH = 2444238.5d;
//    Ecliptic longitude of the Sun at epoch 1980.0
    const ECLIPTIC_LONGITUDE_EPOCH = 278.833540d;
//    Ecliptic longitude of the Sun at perigee
    const ECLIPTIC_LONGITUDE_PERIGEE = 282.596403d;
//    Eccentricity of Earth's orbit
    const ECCENTRICITY = 0.016718d;
//    Elements of the Moon's orbit, epoch 1980.0

//    Moon's mean longitude at the epoch
    const MOON_MEAN_LONGITUDE_EPOCH = 64.975464d;
//    Mean longitude of the perigee at the epoch
    const MOON_MEAN_PERIGEE_EPOCH = 349.383063d;
//    Synodic month (new Moon to new Moon), in days
//var lunardays = 29.53058770576d;
    const SYNODIC_MONTH = 29.53058868d;
    
    const EPSILON = 1e-6;
    
// Precision used when describing the moon's phase in textual format,
// in phase_string().
//	const PRECISION = 0.03;
//	const PRECISION = 0.00000005d;
//	const PRECISION = 0.002;
	const PRECISION = 0.02d;
	const NEW =   0.0 / 4.0d;
	const FIRST = 1.0 / 4.0d;
	const FULL = 2.0 / 4.0d;
	const LAST = 3.0 / 4.0d;
	const NEXTNEW = 4.0 / 4.0d;

	function phaseDescription(south, p) {
   	var phaseStrings = [
       	[NEW + PRECISION, "New"],
       	[FIRST - PRECISION, "Waxing crescent"],
       	[FIRST + PRECISION, "First quarter"],
       	[FULL - PRECISION, "Waxing gibbous"],
       	[FULL + PRECISION, "Full"],
       	[LAST - PRECISION, "Waning gibbous"],
       	[LAST + PRECISION, "Third quarter"],
       	[NEXTNEW - PRECISION, "Waning crescent"],
       	[NEXTNEW + PRECISION, "New"] ];
	var phaseBitmaps = [
		Rez.Drawables.Phase00,
		Rez.Drawables.Phase17,
		Rez.Drawables.Phase26,
		Rez.Drawables.Phase35,
		Rez.Drawables.Phase44,
		Rez.Drawables.Phase53,
		Rez.Drawables.Phase62,
		Rez.Drawables.Phase71
		];

		var i;
		var s = phaseStrings.size();
		for (i = 0; i < s; i++) {
			if (p > phaseStrings[i][0]) {
				continue;
			}
			break;
		}

		var	n = (south ? (s-1-i) : i) % (s-1);
	    return [phaseStrings[i][1], phaseBitmaps[n]];
	}

	function myMod(a, b) {
		var d = (b < 0) ? -b : b;
		var m = (a - ((a / d).toLong() * d));
		var r = ((m < 0) ? (d + m) : m);
		return ((b < 0) ? (r + b) : r);
	}

	function jdn(year, month, day, hour, min, sec) {
		var options = {
			:year   => year,
			:month  => month,
			:day    => day,
			:hour   => hour,
			:minute => min,
			:second => sec
		};
		return Gregorian.moment(options).value() / 86400.0d + 2440587.5d;     // (seconds /(seconds per day)) + julian date of epoch       2440587.5 / 86400 = 28,24753472222 Days
	}

/*  TRUEPHASE  --  Given a K value used to determine the mean phase of
                   the new moon, and a phase selector (0.0, 0.25, 0.5,
                   0.75), obtain the true, corrected phase time.  */

function truephase(k) {
    var t = k / 1236.85d;             /* Time in Julian centuries from
                                       1900 January 0.5 */
    var t2 = t * t;                     /* Square for frequent use */
    var t3 = t2 * t;                    /* Cube for frequent use */
    var pt = 2415020.75933d              /* Mean time of phase */
         + SYNODIC_MONTH * k
         + 0.0001178d * t2
         - 0.000000155d * t3
         + 0.00033d * Math.sin(Math.toRadians(166.56d + 132.87d * t - 0.009173d * t2));

    var m = 359.2242d                    /* Sun's mean anomaly */
        + 29.10535608d * k
        - 0.0000333d * t2
        - 0.00000347d * t3;
    var mprime = 306.0253d               /* Moon's mean anomaly */
        + 385.81691806d * k
        + 0.0107306d * t2
        + 0.00001236d * t3;
    var f = 21.2964d                     /* Moon's argument of latitude */
        + 390.67050646d * k
        - 0.0016528d * t2
        - 0.00000239d * t3;

       /* Corrections for New and Full Moon */

       pt +=     (0.1734d - 0.000393d * t) * Math.sin(Math.toRadians(m))
                + 0.0021d * Math.sin(Math.toRadians(2 * m))
                - 0.4068d * Math.sin(Math.toRadians(mprime))
                + 0.0161d * Math.sin(Math.toRadians(2 * mprime))
                - 0.0004d * Math.sin(Math.toRadians(3 * mprime))
                + 0.0104d * Math.sin(Math.toRadians(2 * f))
                - 0.0051d * Math.sin(Math.toRadians(m + mprime))
                - 0.0074d * Math.sin(Math.toRadians(m - mprime))
                + 0.0004d * Math.sin(Math.toRadians(2 * f + m))
                - 0.0004d * Math.sin(Math.toRadians(2 * f - m))
                - 0.0006d * Math.sin(Math.toRadians(2 * f + mprime))
                + 0.0010d * Math.sin(Math.toRadians(2 * f - mprime))
                + 0.0005d * Math.sin(Math.toRadians(m + 2 * mprime));
    return pt;
}

function getAge(year, month, day, hour, min, sec) {
    var mtime;
    var phaset = new [3];
    
	var now = jdn(year, month, day, hour, min, sec);
    var k1 = Math.floor((year - 1900) * 12.3685).toNumber();// - 4;
    var minx = 0;
    for (var l = 0; true; l++) {
        mtime = truephase(k1);
        if (mtime >= now) {
            minx++;
        }
        phaset[minx] = mtime;
        if (mtime > now) {
            minx++;
            break;
        }
        k1 += 1;
    }
    
    return [now - phaset[0], phaset[1] - phaset[0]]; // age in days, month in days
}

	function fixangle (a) { return (a - 360.0d * Math.floor(a/360.0d)); }

	function kepler(m, ecc) {
/* Solve the equation of Kepler. */
    	m = Math.toRadians(m);
    	var e = m;
    	while (true) {
	        var delta = e - ecc * Math.sin(e) - m;
        	e = e - delta / (1.0d - ecc * Math.cos(e));
	
        	if (delta.abs() <= EPSILON) {
	            break;
        	}
		}
	    return e;
	}

	function getPhase(year, month, day, hour,  min, sec) {
/* Calculate phase of moon as a fraction

    The argument is the time for which the phase is requested,
    expressed as a moment.

    Returns a dictionary containing the terminator phase angle as a
    percentage of a full circle (i.e., 0 to 1), the illuminated
    fraction of the Moon's disc, the Moon's age in days and fraction,
    the distance of the Moon from the centre of the Earth, and the
    angular diameter subtended by the Moon as seen by an observer at
    the centre of the Earth. */

//    date within the epoch
// 	  eg March 29 2021 (Monday) is JDN 2459303.xxx
		var phase_date_jdn = jdn(year, month, day, hour, min, sec);
		var pday = phase_date_jdn - EPOCH;
//    Mean anomaly of the Sun
	    var N = fixangle((360.0d/365.2422d) * pday);
//    Convert from perigee coordinates to epoch 1980
    	var M = fixangle(N + ECLIPTIC_LONGITUDE_EPOCH - ECLIPTIC_LONGITUDE_PERIGEE);

//    Solve Kepler's equation
    	var Ec = kepler(M, ECCENTRICITY);
    	Ec = Math.sqrt((1 + ECCENTRICITY) / (1 - ECCENTRICITY)) * Math.tan(Ec/2.0);
//    True anomaly
    	Ec = 2 * Math.toDegrees(Math.atan(Ec));
//    Suns's geometric ecliptic longuitude
    	var lambda_sun = fixangle(Ec + ECLIPTIC_LONGITUDE_PERIGEE);

//    *********** Calculation of the Moon's position *********
//    Moon's mean longitude
	    var moon_longitude = fixangle(13.1763966d * pday + MOON_MEAN_LONGITUDE_EPOCH);
//    Moon's mean anomaly
	    var MM = fixangle(moon_longitude - 0.1114041d * pday - MOON_MEAN_PERIGEE_EPOCH);
//    Moon's ascending node mean longitude
//    MN = fixangle(NODE_MEAN_LONGITUDE_EPOCH - 0.0529539 * pday);
	    var evection = 1.2739d * Math.sin(Math.toRadians(2*(moon_longitude - lambda_sun) - MM));
//    Annual equation
	    var annual_eq = 0.1858d * Math.sin(Math.toRadians(M));
//    Correction term
	    var A3 = 0.37d * Math.sin(Math.toRadians(M));
    	var MmP = MM + evection - annual_eq - A3;
//    Correction for the equation of the centre
	    var mEc = 6.2886d * Math.sin(Math.toRadians(MmP));
//    Another correction term
	    var A4 = 0.214d * Math.sin(Math.toRadians(2.0d * MmP));
//    Corrected longitude
	    var lP = moon_longitude + evection + mEc - annual_eq + A4;
//    Variation
	    var variation = 0.6583d * Math.sin(Math.toRadians(2*(lP - lambda_sun)));
//    True longitude
	    var lPP = lP + variation;

 // *******    Calculation of the phase of the Moon *******
//    Age of the Moon, in degrees
	    var moon_age = lPP - lambda_sun;
//    Phase of the Moon
	    var moon_phase = fixangle(moon_age) / 360.0d;
	    var moon_illum = (1 - Math.cos(Math.toRadians(moon_age))) / 2.0;
	    var res = {
    	    "phase" =>  moon_phase,
        	"illuminated" => moon_illum,
        	"age"=> getAge(year, month, day, hour,  min, sec)[0]
		};

	    return res;
	}
}
