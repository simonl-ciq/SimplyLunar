using Toybox.Application;

class SimplyLunarApp extends Application.AppBase {

    function initialize() {
        AppBase.initialize();
    }

    // onStart() is called on application start up
    function onStart(state) {
    }

    // onStop() is called when your application is exiting
    function onStop(state) {
    }

    // Return the initial view of your application here
    function getInitialView() {
        var LunarView = new SimplyLunarView();
        return [ LunarView, new SimplyLunarDelegate(LunarView) ];
    }


    // New app settings have been received so trigger a UI update
/*    
    function onSettingsChanged() {
    	if (LunarView != null) {
			LunarView.getDisplayType();
    	    LunarView.requestUpdate();
    	}
    	if (GlanceView != null) {
			GlanceView.getDisplayType();
    	    GlanceView.requestUpdate();
    	}
    }
*/
(:glance)
    function getGlanceView() {
        return [ new SimplyLunarGlanceView() ];
    }

}