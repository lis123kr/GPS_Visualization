from matplotlib.widgets import AxesWidget
from numpy import round
class ZoomSlider(AxesWidget):
    """
    A slider representing a floating point range.
    Create a slider from *valmin* to *valmax* in axes *ax*. For the slider to
    remain responsive you must maintain a reference to it. Call
    :meth:`on_changed` to connect to the slider event.
    Attributes
    ----------
    val : float
        Slider value.
    """
    def __init__(self, ax, label, valmin, valmax, valleft=0.5, valright=1.0, 
                 closedmin=True, closedmax=True, slidermin=None, slidermax=None, 
                 valfmt='%1.2f', dragging=True, valstep=None, **kwargs):
        """
        Parameters
        ----------
        ax : Axes
            The Axes to put the slider in.
        label : str
            Slider label.
        valmin : float
            The minimum value of the slider.
        valmax : float
            The maximum value of the slider.
        closedmin : bool, optional, default: True
            Indicate whether the slider interval is closed on the bottom.
        closedmax : bool, optional, default: True
            Indicate whether the slider interval is closed on the top.
        slidermin : Slider, optional, default: None
            Do not allow the current slider to have a value less than
            the value of the Slider `slidermin`.
        slidermax : Slider, optional, default: None
            Do not allow the current slider to have a value greater than
            the value of the Slider `slidermax`.
        valfmt    : the format string for formatting the slider text
        dragging : bool, optional, default: True
            If True the slider can be dragged by the mouse.
        valstep : float, optional, default: None
            If given, the slider will snap to multiples of `valstep`.
        Notes
        -----
        Additional kwargs are passed on to ``self.poly`` which is the
        :class:`~matplotlib.patches.Rectangle` that draws the slider
        knob.  See the :class:`~matplotlib.patches.Rectangle` documentation for
        valid property names (e.g., `facecolor`, `edgecolor`, `alpha`).
        """
        AxesWidget.__init__(self, ax)

        if slidermin is not None and not hasattr(slidermin, 'val'):
            raise ValueError("Argument slidermin ({}) has no 'val'"
                             .format(type(slidermin)))
        if slidermax is not None and not hasattr(slidermax, 'val'):
            raise ValueError("Argument slidermax ({}) has no 'val'"
                             .format(type(slidermax)))
        self.closedmin = closedmin
        self.closedmax = closedmax
        self.slidermin = slidermin
        self.slidermax = slidermax
        self.drag_active = False
        self.valmin = valmin
        self.valmax = valmax
        self.valstep = valstep
        self.valfmt = valfmt
        valleft = self._value_in_bounds(valleft)
        valright = self._value_in_bounds(valright)
        self.valleft = valleft
        self.valright = valright

    
        self.poly = ax.axvspan(valleft, valright, 0, 1, **kwargs)
        # self.polyright = ax.axvspan(valright, valmax, 0, 1, **kwargs)

        ax.set_yticks([])
        ax.set_xlim((valmin, valmax))
        ax.set_xticks([])
        ax.set_navigate(False)

        self.connect_event('button_press_event', self._update)
        self.connect_event('button_release_event', self._update)
        if dragging:
            self.connect_event('motion_notify_event', self._update)
        self.label = ax.text(-0.02, 0.5, label, transform=ax.transAxes,
                             verticalalignment='center',
                             horizontalalignment='right')
        self.lefttext = ax.text(0.0, -0.4, valfmt % valleft,
                               transform=ax.transAxes,
                               verticalalignment='center',
                               horizontalalignment='left')
        self.righttext = ax.text(0.95, -0.4, valfmt % valright,
                               transform=ax.transAxes,
                               verticalalignment='center',
                               horizontalalignment='left')

        self.cnt = 0
        self.observers = {}

        self.set_leftval(valleft)
        self.set_rightval(valright)

    def _value_in_bounds(self, val):
        """ Makes sure self.val is with given bounds."""
        if self.valstep:
            val = round((val - self.valmin)/self.valstep)*self.valstep
            val += self.valmin

        if val <= self.valmin:
            if not self.closedmin:
                return
            val = self.valmin
        elif val >= self.valmax:
            if not self.closedmax:
                return
            val = self.valmax

        if self.slidermin is not None and val <= self.slidermin.val:
            if not self.closedmin:
                return
            val = self.slidermin.val

        if self.slidermax is not None and val >= self.slidermax.val:
            if not self.closedmax:
                return
            val = self.slidermax.val
        return val

    def _update(self, event):
        """update the slider position"""
        if self.ignore(event):
            return

        if event.button != 1:
            return

        if event.name == 'button_press_event' and event.inaxes == self.ax:
            self.drag_active = True
            event.canvas.grab_mouse(self.ax)

        if not self.drag_active:
            return

        elif ((event.name == 'button_release_event') or
              (event.name == 'button_press_event' and
               event.inaxes != self.ax)):
            self.drag_active = False
            event.canvas.release_mouse(self.ax)
            return
        val = self._value_in_bounds(event.xdata)
        if (val is not None) and (val != self.valleft) and abs(val - self.valleft) < abs(val - self.valright):
            self.set_leftval(val)
        elif (val is not None) and (val != self.valright) and abs(val - self.valleft) > abs(val - self.valright):
            self.set_rightval(val)

    def set_leftval(self, val):
        """
        Set slider value to *val*
        Parameters
        ----------
        val : float
        """
        xy = self.poly.xy
        xy[0] = val, 0
        xy[1] = val, 1
        self.poly.xy = xy
        self.lefttext.set_text(self.valfmt %val)
        if self.drawon:
            self.ax.figure.canvas.draw_idle()
        self.valleft = val
        if not self.eventson:
            return
        for cid, func in self.observers.items():
            func(val)

    def set_rightval(self, val):
        """
        Set slider value to *val*
        Parameters
        ----------
        val : float
        """
        xy = self.poly.xy
        xy[2] = val, 1
        xy[3] = val, 0
        self.poly.xy = xy
        self.righttext.set_text(self.valfmt %val)
        if self.drawon:
            self.ax.figure.canvas.draw_idle()
        self.valright = val
        if not self.eventson:
            return
        for cid, func in self.observers.items():
            func(val)

    def on_changed(self, func):
        """
        When the slider value is changed call *func* with the new
        slider value
        Parameters
        ----------
        func : callable
            Function to call when slider is changed.
            The function must accept a single float as its arguments.
        Returns
        -------
        cid : int
            Connection id (which can be used to disconnect *func*)
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        """
        Remove the observer with connection id *cid*
        Parameters
        ----------
        cid : int
            Connection id of the observer to be removed
        """
        try:
            del self.observers[cid]
        except KeyError:
            pass

    def reset(self):
        """Reset the slider to the initial value"""
        if (self.valleft != self.valmin):
            self.set_leftval(self.valmin)
        if (self.valright != self.valmax):
            self.set_rightval(self.valmax)