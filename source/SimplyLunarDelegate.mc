using Toybox.WatchUi as Ui;
using Toybox.System as Sys;

class SimplyLunarDelegate extends Ui.BehaviorDelegate {
    /* Initialize and get a reference to the view, so that
     * user iterations can call methods in the main view. */
	var SLView;
     
    function initialize(view) {
        Ui.BehaviorDelegate.initialize();
        SLView = view;
    }

    function onSelect() {
        Ui.requestUpdate();
        return true;
    }

    function onMenu() {
        SLView.setPage((SLView.Page + 1) % 2);
        Ui.requestUpdate();
        return true;
    }
    
}
