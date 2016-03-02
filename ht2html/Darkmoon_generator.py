"""
Generates the ht2html documentation style
"""
import os
import time

from Skeleton import Skeleton
from Sidebar import Sidebar, BLANKCELL
from Banner import Banner
from HTParser import HTParser
from LinkFixer import LinkFixer

sitelinks=[]
# You can add entries via file banner.h
# if you want
if os.path.exists('banner.h'):
    banneritems=file('banner.h').readlines()
    for item in banneritems:
        link, name = item.split("\t")
        sitelinks.append((link, name))

class Darkmoon_generator(Skeleton, Sidebar, Banner):
    AUTHOR = 'darkmoon'
    EMAIL = 'darkmoon@altervista.org'

    def __init__(self, file, rootdir, relthis):
        root, ext = os.path.splitext(file)
        html = root + '.html'
        p = self.__parser = HTParser(file, self.AUTHOR, self.EMAIL)
        f = self.__linkfixer = LinkFixer(html, rootdir, relthis)
        self.__body = None
        self.__cont = None
        # Calculate the sidebar links, adding a few of our own.
        self.__d = {'rootdir': rootdir}
        p.process_sidebar()
        p.sidebar.append(BLANKCELL)
        # It is important not to have newlines between the img tag and the end
        # end center tags, otherwise layout gets messed up.
        p.sidebar.append(('http://www.python.org/', '''\
<img class="side" alt="[Python Powered]" border="0"
src="img/PythonPoweredSmall.png">''' % self.__d))
        self.__linkfixer.massage(p.sidebar, self.__d)
        Sidebar.__init__(self, p.sidebar)
        copyright = self.__parser.get('copyright', '1999-%d' %
                                      time.localtime()[0])
        p.sidebar.append((None, '&copy; ' + copyright))
        p.sidebar.append(('http://www.python.org/psf/',
                          'Python Software Foundation'))
        p.sidebar.append((None, '<hr>'))
        p.sidebar.append((None, '''<script type="text/javascript" \
        src="http://www.altervista.org/js_tags/contatore.js">
        </script>
        '''))
        p.sidebar.append((None, '''last update'''))
        p.sidebar.append((None, time.asctime()))
        
        # Fix up our site links, no relthis because the site links are
        # relative to the root of our web pages.
        sitelink_fixer = LinkFixer(f.myurl(), rootdir)
        sitelink_fixer.massage(sitelinks, self.__d, aboves=1)
        Banner.__init__(self, sitelinks)
        # kludge!
        for i in range(len(p.sidebar)-1, -1, -1):
            if p.sidebar[i] == 'Email Us':
                p.sidebar[i] = 'Email me'
                break

    def get_title(self):
        return self.__parser.get('title')

    def get_sidebar(self):
        if self.__parser.get('wide-page', 'no').lower() == 'yes':
            return None
        return Sidebar.get_sidebar(self)

    def get_banner(self):
        return Banner.get_banner(self)

    def get_banner_attributes(self):
        return 'cellspacing="0" cellpadding="0"'

    def get_corner(self):
        # It is important not to have newlines between the img tag and the end
        # anchor and end center tags, otherwise layout gets messed up
        return '''<a class="corner" href="index.html"><img class="corner"
        alt='-=[ `/\- DARKsiteoftheMOON `/\- ]=-' border="0"
         src="img/corner.png"></a>''' % self.__d

    def get_corner_bgcolor(self):
        # this may not be 100% correct.  it uses PIL to get the RGB values at
        # the corners of the image and then takes a vote as to the most likely
        # value.  Some images may be `bizarre'.  See .../pics/backgrounds.py
        return '#4fa445'

    def get_body(self):
        self.__grokbody()
        return self.__body

    def get_cont(self):
        self.__grokbody()
        return self.__cont

    def __grokbody(self):
        if self.__body is None:
            text = self.__parser.fp.read()
            i = text.find('<!--table-stop-->')
            if i >= 0:
                self.__body = text[:i]
                self.__cont = text[i+17:]
            else:
                # there is no wide body
                self.__body = text
