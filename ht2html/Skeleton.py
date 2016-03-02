"""Skeleton class.

Should be sub-classed to provide basic generation of able-contained HTML
document.
"""

from ht2html import __version__

import os
import sys
import time
from cStringIO import StringIO


class Skeleton:
    #
    # for sub-classes to override
    #

    def get_banner(self):
        """Returns HTML for the top banner, or None if no banner.
        """
        return None

    def get_left_sidebar(self):
        """Returns HTML for the left sidebar or None.
        """
        return None
    # for backwards compatibility
    get_sidebar = get_left_sidebar

    def get_right_sidebar(self):
        """Returns HTML for the right sidebar or None.
        """
        return None

    def get_banner_width(self):
        """HTML `width' of banner column as a percentage.

        Should be a string that does not include the percent sign (e.g. "90").
        This affects the column containing both the banner and body text (if
        they exist).
        """
        return '90'

    def get_corner(self):
        """Returns HTML for the upper-left corner or None.

        Note that if get_banner() and get_sidebar() both return None, then
        get_corner() is ignored.  Also if both get_banner() and get_sidebar()
        return a string, but get_corner() returns None, the smallest blank
        corner possible is emitted.
        """
        return None

    def get_body(self):
        """Returns HTML for the internal document body.

        Note that using this method causes the get_sidebar() -- if there is
        one -- to run the full height of the page.  If you don't want this,
        then make sure get_cont() returns a string.
        """
        return '<b>Intentionally left blank</b>'

    def get_cont(self):
        """Returns HTML for the continuation of the body.

        When this method returns a string, and there is a get_sidebar(), the
        continuation text appears below the get_sidebar() and get_body() at
        the full width of the screen.  If there is no get_sidebar(), then
        having a get_cont() is pointless.
        """
        return None

    def get_title(self):
        """Return the title of the page.  Required."""
        return 'Intentionally left blank'

    def get_meta(self):
        """Return extra meta-data.  Must be a string."""
        import __main__
        return '<meta name="generator" content="HT2HTML/%s">' \
               % __main__.__version__

    def get_headers(self):
        """Return extra header information.  Must be a string."""
        return ''

    def get_corner_bgcolor(self):
        """Return the background color for the corner"""
        return self.get_lightshade()

    def get_body_attributes(self):
        """Return extra attributes for the body start tag."""
        # These are not supported in HTML, but are needed for
        # Netscape 4
        return 'marginwidth="0" marginheight="0"'

    def get_banner_attributes(self):
        """Return extra attributes for the TABLE in the banner."""
        return 'cellspacing="0" cellpadding="2"'

    def get_charset(self):
        """Return charset of pages"""
        return 'it-ascii'

    # Style sheets
    def get_stylesheet(self):
        """Return filename of CSS stylesheet."""
        return 'css/style.css'

    def get_stylesheet_pi(self):
        s = self.get_stylesheet()
        if s:
            return '<?xml-stylesheet href="%s" type="%s"?>\n' \
                   % (s, self.get_stylesheet_type(s))
        else:
            return ''

    def get_stylesheet_type(self, filename):
        ext = os.path.splitext(filename)[1]
        if ext == ".css":
            return "text/css"
        elif ext in (".xsl", ".xslt"):
            return "text/xslt"
        else:
            raise ValueError("unknown stylesheet language")

    def get_style(self):
        """Return the style sheet for this document"""
        s = self.body_style()
        if s:
            return 'body { %s }' % self.body_style()
        else:
            return ''

    def body_style(self):
        if self.get_stylesheet():
            # If there's an external stylesheet, rely on that for the body.
            return ''
        else:
            return 'margin: 0px;'

    # Call this method
    def makepage(self):
        banner = self.get_banner()
        sidebar = self.get_sidebar()
        corner = self.get_corner()
        body = self.get_body()
        cont = self.get_cont()
        html = StringIO()
        stdout = sys.stdout
        closed = 0
        try:
            sys.stdout = html
            self.__do_head()
            self.__start_body()
            print '<!-- start of page table -->'
            print ('<table summary="Page" width="100%" border="0"'
                   ' cellspacing="0" cellpadding="0">')
            if banner is not None:
                print '<!-- start of banner row -->'
                print '<tr>'
                if corner is not None:
                    self.__do_corner(corner)
                print '<!-- start of banner -->'
                print '<td width="%s%%" class="bannerdark">' %\
                      self.get_banner_width()
		print banner
                print '</td><!-- end of banner -->'
                print '</tr><!-- end of banner row -->'
            # if there is a body but no sidebar, then we'll just close the
            # table right here and put the body (and any cont) in the full
            # page.  if there is a sidebar but no body, then we still create
            # the new row and just populate the body cell with a non-breaking
            # space.  Watch out though because we don't want to close the
            # table twice
            if sidebar is None:
                print '</table><!-- end of page table -->'
                closed = 1
            else:
                print '<tr><!-- start of sidebar/body row -->'
                self.__do_sidebar(sidebar)
            if body is not None:
                if closed:
                    print body
                else:
                    self.__do_body(body)
            if not closed:
                print '</tr><!-- end of sidebar/body row -->'
                print '</table><!-- end of page table -->'
            if cont is not None:
                self.__do_cont(cont)
            self.__finish_all()
        finally:
            sys.stdout = stdout
        return html.getvalue()

    def __do_corner(self, corner):
        print '<!-- start of corner cells -->'
        print '<td width="150" valign="middle" class="bannerdark">'
        # it is important not to have a newline between the corner text and
        # the table close tag, otherwise layout is messed up
        if corner is None:
            print '&nbsp;',
        else:
            print corner,
        print '</td>'
        print '<!-- end of corner cells -->'

    def __do_sidebar(self, sidebar):
        print '<!-- start of sidebar cells -->'
        print '<td class="sidelight" width="150" valign="top" >'
        print sidebar
        print '</td>'
#
#  Space unneeded (use stylesheet instead ...)
#        print '<td class="body" width="1">&nbsp;</td><!--spacer-->'
        print '<!-- end of sidebar cell -->'

    def __do_body(self, body):
        print '<!-- start of body cell -->'
        print '<td valign="top" width="%s%%" class="body"><br>' % (
            self.get_banner_width())
        print """<div id="container">
<b class="rtop">
<b class="r1"></b> <b class="r2"></b> <b class="r3"></b> <b class="r4"></b>
</b>
"""
        print body
        print """<b class="rbottom">
   <b class="r4"></b> <b class="r3"></b> <b class="r2"></b> <b class="r1"></b>
 </b>
 </div>
        """
        print '<p>&nbsp;</p></td><!-- end of body cell -->'

    def __do_cont(self, cont):
        print '<div class="body">'
        print '<div class="continuation">'
        print '<!-- start of continued wide-body text -->'
        print cont
        print '<!-- end of continued wide-body text -->'
        print '</div>'
        print '</div>'

    def __do_head(self):
        """Return the HTML <head> stuff."""
        print '''\
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
%(xmlstyle)s<html>
<!-- THIS PAGE IS AUTOMATICALLY GENERATED.  DO NOT EDIT. -->
<!-- %(time)s -->
<!-- USING HT2HTML %(version)s -->
<!-- SEE http://ht2html.sf.net -->
<!-- User-specified headers:
Title: %(title)s
%(headers)s
-->

<head>
<title>%(title)s</title>
<meta http-equiv="Content-Type" content="text/html; charset=%(charset)s">
%(meta)s
%(style)s
</head>''' % {'title'   : self.get_title(),
              'headers' : self.get_headers(),
              'meta'    : self.get_meta(),
              'time'    : time.ctime(time.time()),
              'version' : __version__,
              'charset' : self.get_charset(),
              'style'   : self.__do_styles(),
              'xmlstyle': self.get_stylesheet_pi(),
              }

    def __do_styles(self):
        # assemble all the style information we have to produce the
        # appropriate LINK and STYLE elements
        stylesheet = self.get_stylesheet()
        localstyle = self.get_style()
        s = ''
        if stylesheet and stylesheet.strip():
            stylesheet = stylesheet.strip()
            type = self.get_stylesheet_type(stylesheet)
            s = '<link rel="STYLESHEET" href="%s" type="%s">' \
                % (stylesheet, type)
        if localstyle and localstyle.strip():
            localstyle = '<style type="text/css">\n%s\n</style>' \
                         % localstyle.strip()
            if stylesheet:
                s = s + "\n" + localstyle
            else:
                s = localstyle
        return s

    def __start_body(self):
        print '''\
<body class=body>'''

    def __finish_all(self):
        print '<p class="time">Last update: %s</p>' % time.asctime()
        print '''<script src="http://www.google-analytics.com/urchin.js" type="text/javascript">
</script>
<script type="text/javascript">
_uacct = "UA-2587939-5";
urchinTracker();
</script>'''
        print '''<center><script type="text/javascript"><!--
google_ad_client = "pub-0988154404168132";
google_ad_width = 728;
google_ad_height = 90;
google_ad_format = "728x90_as";
google_ad_type = "text_image";
google_ad_channel = "";
google_color_border = "000000";
google_color_bg = "000000";
google_color_link = "FFFFFF";
google_color_text = "CCCCCC";
google_color_url = "999999";
google_ui_features = "rc:0";
//-->
</script>
<script type="text/javascript"
  src="http://pagead2.googlesyndication.com/pagead/show_ads.js">
</script></center>'''
        print '</body></html>'

# test script
class _Skeleton(Skeleton):
    def get_banner(self):
        return '<b>The Banner</b>'

    def get_sidebar(self):
        return '''<ul class=side><li class=side>Sidebar line 1
        <li class=side>Sidebar line 2
        <li class=side>Sidebar line 3
        </ul>'''

    def get_corner(self):
        return '<center><em class=side>CORNER</em></center>'

    def get_body(self):
        return 'intentionally left blank ' * 10

    def get_cont(self):
        return 'wide stuff ' * 10

    def get_corner_bgcolor(self):
        return 'yellow'

    def get_banner_width(self):
        return "80"


if __name__ == '__main__':
    t = _Skeleton()
    print t.makepage()
