# This file is generated by wxPython's SIP generator.  Do not edit by hand.
# 
# Copyright: (c) 2013 by Total Control Software
# License:   wxWindows License

from ._adv import *

import wx

EVT_DATE_CHANGED = wx.PyEventBinder( wxEVT_DATE_CHANGED, 1 )
EVT_TIME_CHANGED = wx.PyEventBinder( wxEVT_TIME_CHANGED, 1 )

DateEvent.PyGetDate = wx.deprecated(DateEvent.GetDate, 'Use GetDate instead.')
DateEvent.PySetDate = wx.deprecated(DateEvent.SetDate, 'Use SetDate instead.')

EVT_CALENDAR =                 wx.PyEventBinder( wxEVT_CALENDAR_DOUBLECLICKED, 1)
EVT_CALENDAR_SEL_CHANGED =     wx.PyEventBinder( wxEVT_CALENDAR_SEL_CHANGED, 1)
EVT_CALENDAR_WEEKDAY_CLICKED = wx.PyEventBinder( wxEVT_CALENDAR_WEEKDAY_CLICKED, 1)
EVT_CALENDAR_PAGE_CHANGED =    wx.PyEventBinder( wxEVT_CALENDAR_PAGE_CHANGED, 1)

EVT_CALENDAR_DAY =             wx.PyEventBinder( wxEVT_CALENDAR_DAY_CHANGED, 1)
EVT_CALENDAR_MONTH =           wx.PyEventBinder( wxEVT_CALENDAR_MONTH_CHANGED, 1)
EVT_CALENDAR_YEAR =            wx.PyEventBinder( wxEVT_CALENDAR_YEAR_CHANGED, 1)

CalendarCtrl.PyGetDate = wx.deprecated(CalendarCtrl.GetDate, 'Use GetDate instead.')
CalendarCtrl.PySetDate = wx.deprecated(CalendarCtrl.SetDate, 'Use SetDate instead.')
CalendarCtrl.PySetDateRange = wx.deprecated(CalendarCtrl.SetDateRange, 'Use SetDateRange instead.')

GenericCalendarCtrl.PyGetDate = wx.deprecated(GenericCalendarCtrl.GetDate, 'Use GetDate instead.')
GenericCalendarCtrl.PySetDate = wx.deprecated(GenericCalendarCtrl.SetDate, 'Use SetDate instead.')
GenericCalendarCtrl.PySetDateRange = wx.deprecated(GenericCalendarCtrl.SetDateRange, 'Use SetDateRange instead.')

EVT_HYPERLINK = wx.PyEventBinder( wxEVT_HYPERLINK, 1 )

# deprecated wxEVT alias
wxEVT_COMMAND_HYPERLINK  = wxEVT_HYPERLINK

EVT_TASKBAR_MOVE = wx.PyEventBinder (         wxEVT_TASKBAR_MOVE )
EVT_TASKBAR_LEFT_DOWN = wx.PyEventBinder (    wxEVT_TASKBAR_LEFT_DOWN )
EVT_TASKBAR_LEFT_UP = wx.PyEventBinder (      wxEVT_TASKBAR_LEFT_UP )
EVT_TASKBAR_RIGHT_DOWN = wx.PyEventBinder (   wxEVT_TASKBAR_RIGHT_DOWN )
EVT_TASKBAR_RIGHT_UP = wx.PyEventBinder (     wxEVT_TASKBAR_RIGHT_UP )
EVT_TASKBAR_LEFT_DCLICK = wx.PyEventBinder (  wxEVT_TASKBAR_LEFT_DCLICK )
EVT_TASKBAR_RIGHT_DCLICK = wx.PyEventBinder ( wxEVT_TASKBAR_RIGHT_DCLICK )
EVT_TASKBAR_CLICK =  wx.PyEventBinder (       wxEVT_TASKBAR_CLICK )
EVT_TASKBAR_BALLOON_TIMEOUT = wx.PyEventBinder ( wxEVT_TASKBAR_BALLOON_TIMEOUT )
EVT_TASKBAR_BALLOON_CLICK = wx.PyEventBinder ( wxEVT_TASKBAR_BALLOON_CLICK )

SPLASH_CENTER_ON_PARENT = SPLASH_CENTRE_ON_PARENT
SPLASH_CENTER_ON_SCREEN = SPLASH_CENTRE_ON_SCREEN
SPLASH_NO_CENTER = SPLASH_NO_CENTRE

EVT_SASH_DRAGGED = wx.PyEventBinder( wxEVT_SASH_DRAGGED, 1 )
EVT_SASH_DRAGGED_RANGE = wx.PyEventBinder( wxEVT_SASH_DRAGGED, 2 )

EVT_QUERY_LAYOUT_INFO = wx.PyEventBinder( wxEVT_QUERY_LAYOUT_INFO )
EVT_CALCULATE_LAYOUT = wx.PyEventBinder( wxEVT_CALCULATE_LAYOUT )

EVT_WIZARD_BEFORE_PAGE_CHANGED  = wx.PyEventBinder( wxEVT_WIZARD_BEFORE_PAGE_CHANGED, 1)
EVT_WIZARD_PAGE_CHANGED  = wx.PyEventBinder( wxEVT_WIZARD_PAGE_CHANGED, 1)
EVT_WIZARD_PAGE_CHANGING = wx.PyEventBinder( wxEVT_WIZARD_PAGE_CHANGING, 1)
EVT_WIZARD_CANCEL        = wx.PyEventBinder( wxEVT_WIZARD_CANCEL, 1)
EVT_WIZARD_HELP          = wx.PyEventBinder( wxEVT_WIZARD_HELP, 1)
EVT_WIZARD_FINISHED      = wx.PyEventBinder( wxEVT_WIZARD_FINISHED, 1)
EVT_WIZARD_PAGE_SHOWN    = wx.PyEventBinder( wxEVT_WIZARD_PAGE_SHOWN, 1)

