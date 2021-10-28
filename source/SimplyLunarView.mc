using Toybox.WatchUi as Ui;
using Toybox.Graphics as Gfx;
using Toybox.Application as App;
using Toybox.Application.Properties as Props;
using Toybox.Application.Storage as Storage;
using Toybox.Time as Time;
using Toybox.Position as Position;
using Toybox.Time.Gregorian;
using Toybox.System as Sys;

using Astro  as astro;
using Phase  as phase;
/*
using cPhase as cphase;
using Moon as moon;
using AIT as ait;
using TRUE as True;
using IDL as idl;
*/

class SimplyLunarView extends Ui.View {
	var myInfo = null;
	var needGPS = true;

	var Page = 0;

	var Places = 0;
	var Precision = "%.1f";

	var Title = "Lunar Angles";
	var Azimuth = "Azimuth";
	var Elevation = "Elevation";
	var NoGPS = "acquiring GPS";

	var phaseTitle = "Lunar Phase";
	var lunarBitmaps;

    function initialize() {

        View.initialize();
	}

    function onPosition(info) {
    	myInfo = info;
        Ui.requestUpdate();
    }

/*
    function setPrecision() {
		Places = (Places + 1) % 2;
        Storage.setValue("places", Places);
		Precision = "%." + Places.toString() + "f";
	}

    function getPrecision() {
    	var temp;
    	temp = Storage.getValue("places");
        if (temp == null) {
	        Storage.setValue("places", Places);
   		} else {
   			Places = temp;
   		}
		Precision = "%." + Places.toString() + "f";
	}
*/
	
	function setPage (newPage) {
		Page = newPage;
        Storage.setValue("page", Page);
	}

    function getPage() {
    	var temp;
    	temp = Storage.getValue("page");
        if (temp == null) {
	        Storage.setValue("page", Page);
   		} else {
   			Page = temp;
   		}
	}

    // Load your resources here
    function onLayout(dc) {

        setLayout(Rez.Layouts.MainLayout(dc));

		var tmpView = View.findDrawableById("title");
		tmpView.setText(Title);

		var azitView = View.findDrawableById("azititle");
		azitView.setText(Azimuth);

		var eletView = View.findDrawableById("eletitle");
		eletView.setText(Elevation);

    }

    // Called when this View is brought to the foreground. Restore
    // the state of this View and prepare it to be shown. This includes
    // loading resources into memory.
    function onShow() {
		myInfo = Position.getInfo();
        if (myInfo == null || myInfo.accuracy < Position.QUALITY_POOR) {
            Position.enableLocationEvents(Position.LOCATION_CONTINUOUS, method(:onPosition));
		}

		getPage();
/*		
		if (Page == 0) {
			getPrecision();
		}
*/
    }

    // Update the view
    function onUpdate(dc) {

		if (needGPS) {
	    	if (myInfo == null || myInfo.accuracy == null || myInfo.accuracy < Position.QUALITY_POOR) {
		    	myInfo = Position.getInfo();
		    }
			if (myInfo.accuracy != null && myInfo.accuracy != Position.QUALITY_NOT_AVAILABLE && myInfo.position != null) {
				if (myInfo.accuracy >= Position.QUALITY_POOR) {
		            Position.enableLocationEvents(Position.LOCATION_DISABLE, method(:onPosition));
					needGPS = false;
	    		}
	    	}
		}

        if (Page == 0) {
			updateAngle(myInfo, dc);
        // Call the parent onUpdate function to redraw the layout
	        View.onUpdate(dc);
		} else {
		    updatePhase(myInfo, dc);
		}
	}

    function updateAngle(myInfo, dc) {
    	var az = "";
    	var el = "";
    	var azit = Azimuth;
    	var elet = Elevation;

		var azView = View.findDrawableById("azimuth");
		var azitView = View.findDrawableById("azititle");

		var elView = View.findDrawableById("elevation");
		var eletView = View.findDrawableById("eletitle");
		
		var myAccuracy = (!(myInfo has :accuracy) || myInfo.accuracy == null) ? Position.QUALITY_GOOD : myInfo.accuracy;
		if (myAccuracy > Position.QUALITY_NOT_AVAILABLE && myInfo.position != null) {
			var loc = myInfo.position.toDegrees();
			var alt = myInfo.altitude / 1000;
   			var today = Gregorian.utcInfo(Time.now(), Time.FORMAT_SHORT);
			var azimelev = astro.LunarAzEl(today.year, today.month, today.day, today.hour, today.min, today.sec, loc[0], loc[1], (alt < 0) ? 0 : alt);
//	azel = astro.LunarAzEl(2021, 03, 23, 13, 38, 0, 53.32257, -2.6454, -15);

			az = azimelev[0].format(Precision);
			if (az.equals("360")) { az = "0"; }
			else if (az.equals("360.0")) { az = "0.0"; }
			el = azimelev[1].format(Precision);
			if (myAccuracy == Position.QUALITY_LAST_KNOWN) {
				azView.setColor(Gfx.COLOR_LT_GRAY);
				elView.setColor(Gfx.COLOR_LT_GRAY);
			} else {
				azView.setColor(Gfx.COLOR_WHITE);
				elView.setColor(Gfx.COLOR_WHITE);
			}
		} else {
    		azit = NoGPS + " ...";
    		elet = "";
		}

		azitView.setText(azit);
		azView.setText(az);
		eletView.setText(elet);
		elView.setText(el);
    }

/*
hidden var days = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
function getDays(year, month) {
	var len = days[month];
	if (month == 2) {
		if (((year % 400) == 0) || (((year % 4) == 0) && ((year % 100) != 0))) {
			len++;
		}
	}
	return len;
}
*/
(:vivo3)
	function getTitleY() {
        return 12;
	}

(:approachs60)
	function getTitleY() {
        return 15;
	}

(:tinyFont)
	function getTitleY() {
        return 10;
	}

(:mediumFont)
	function getTitleY() {
        return 0;
	}

(:vivo3)
	function getPhaseFonts() {
        return [Graphics.FONT_SMALL, Graphics.FONT_SMALL];
	}

(:approachs60)
	function getPhaseFonts() {
        return [Graphics.FONT_TINY, Graphics.FONT_SMALL];
	}

(:tinyFont)
	function getPhaseFonts() {
        return [Graphics.FONT_TINY, Graphics.FONT_MEDIUM];
	}

(:mediumFont)
	function getPhaseFonts() {
        return [Graphics.FONT_MEDIUM, Graphics.FONT_MEDIUM];
	}

	function updatePhase(myInfo, dc) {
        var hCentre = dc.getWidth() / 2;
   	    var vCentre = dc.getHeight() / 2;

		dc.setColor(Graphics.COLOR_BLACK,Graphics.COLOR_BLACK);
		dc.clear();
		dc.setColor(Graphics.COLOR_WHITE,Graphics.COLOR_TRANSPARENT);

		var fnts = getPhaseFonts();
// Title
		var ty = getTitleY();
        dc.drawText(hCentre, ty, fnts[0], phaseTitle, Graphics.TEXT_JUSTIFY_CENTER);
		
		var today = Gregorian.utcInfo(Time.now(), Time.FORMAT_SHORT);
		var p = phase.getPhase(today.year, today.month, today.day, today.hour, today.min, today.sec);
		var south = (myInfo.accuracy != null &&
			myInfo.accuracy > Position.QUALITY_NOT_AVAILABLE &&
			myInfo.position != null &&
			myInfo.position.toDegrees()[0] < 0 );

		var desc = phase.phaseDescription(south, p["phase"]);
		var illum = (p["illuminated"]*100).format("%.1f") + "%";
		if (illum.equals("100.0%")) { illum = "100%"; }
		var days = (p["age"].format("%.1f"));
		var age = days + " days";// + (days.equals("1") ? "" : "s");

        var bitmap = Ui.loadResource(desc[1]);
		var bw = bitmap.getWidth();
        var bh = bitmap.getHeight();

// Moon
		var bx = hCentre - (bw * 12 / 10);
        var by = vCentre - (bh / 2);
       	dc.drawBitmap(bx, by, bitmap);

// Description
		var dy = vCentre - (bh * 5 / 9) - dc.getFontHeight(Graphics.FONT_MEDIUM);
   		dc.drawText(hCentre, dy, fnts[1], desc[0], Graphics.TEXT_JUSTIFY_CENTER);

// Age
        var ay = vCentre + (bh * 5 / 9);
   	    dc.drawText(hCentre, ay, Graphics.FONT_MEDIUM, age, Graphics.TEXT_JUSTIFY_CENTER);

// Illumination
		var fieldWidth = dc.getTextWidthInPixels("100%", Graphics.FONT_NUMBER_MILD);
		var ix = hCentre + (bw * 12 / 10);// + fieldWidth;// + (bw / 10);
   	    var iy = vCentre - dc.getFontHeight(Graphics.FONT_NUMBER_MILD) / 2;
       	dc.drawText(ix, iy, Graphics.FONT_NUMBER_MILD, illum, Graphics.TEXT_JUSTIFY_RIGHT);
		
	}

    // Called when this View is removed from the screen. Save the
    // state of this View here. This includes freeing resources from
    // memory.
    function onHide() {
        Position.enableLocationEvents(Position.LOCATION_DISABLE, method(:onPosition));
    }

}

(:glance)
class SimplyLunarGlanceView extends Ui.GlanceView {
	var vcentre = 30;
	var hwidth = 191;
	var Page = 0;

	function initialize() {
		GlanceView.initialize();
	}
	
	function onLayout(dc) {
//		vcentre = dc.getFontHeight(Gfx.FONT_SMALL) - 2;
		vcentre = 30;
		hwidth = dc.getWidth();
	}

    function getPage() {
    	var temp;
    	temp = Storage.getValue("page");
        if (temp == null) {
	        Storage.setValue("page", Page);
   		} else {
   			Page = temp;
   		}
	}

	function onUpdate(dc) {
		var az = "";
		var el = "";
	    var	myInfo = Position.getInfo();
//	    	    myInfo.position = Position.parse("53.825564, -2.421976", Position.GEO_DEG);
//	    	    myInfo.position = Position.parse("34.0522, -118.2437", Position.GEO_DEG);
// Tokyo 35.6762, 139.6503

		var today = Gregorian.utcInfo(Time.now(), Time.FORMAT_SHORT);
		getPage();
		if (Page == 0) {
			var myAccuracy = (!(myInfo has :accuracy) || myInfo.accuracy == null) ? Position.QUALITY_GOOD : myInfo.accuracy;
			if (myAccuracy > Position.QUALITY_NOT_AVAILABLE && myInfo.position != null) {
				var loc = myInfo.position.toDegrees();
				var alt = myInfo.altitude / 1000;
				var azimelev = astro.LunarAzEl(today.year, today.month, today.day, today.hour, today.min, today.sec, loc[0], loc[1], (alt < 0) ? 0 : alt);
				az = azimelev[0].format("%.0f");
				if (az.equals("360")) {az = "0";}
				az = "a: "+az+"°";
				el = "e: "+azimelev[1].format("%.0f")+"°";
			}
		} else {
			var p = phase.getPhase(today.year, today.month, today.day, today.hour, today.min, today.sec);
			var illum = (p["illuminated"]*100).format("%.0f");
			var age = p["age"].format("%.0f");
			el = "age: "+ age;
			az = "vis: "+ illum;
		}
		dc.setColor(Graphics.COLOR_BLACK,Graphics.COLOR_BLACK);
		dc.clear();
		dc.setColor(Graphics.COLOR_WHITE,Graphics.COLOR_TRANSPARENT);
        var AppName = Ui.loadResource( Rez.Strings.AppName );
        if (AppName == null ) {AppName = "Simply Lunar";}
		dc.drawText(0, 0, Gfx.FONT_SMALL, AppName, Gfx.TEXT_JUSTIFY_LEFT);
		dc.setColor(Gfx.COLOR_BLUE, Gfx.COLOR_TRANSPARENT);
		if (myInfo.accuracy == Position.QUALITY_LAST_KNOWN) {
			dc.setColor(Gfx.COLOR_LT_GRAY, Gfx.COLOR_TRANSPARENT);
		}
		dc.drawText(0, vcentre, Gfx.FONT_TINY, az, Gfx.TEXT_JUSTIFY_LEFT);
		var x = dc.getTextDimensions(AppName, Gfx.FONT_SMALL)[0];
		dc.drawText(x, vcentre, Gfx.FONT_TINY, el, Gfx.TEXT_JUSTIFY_RIGHT);
	}

}