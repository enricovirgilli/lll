"""Generates the ht2html documentation style
"""
import os
import time

from Skeleton import Skeleton
from Sidebar import Sidebar, BLANKCELL
from Banner import Banner
from HTParser import HTParser
from LinkFixer import LinkFixer



def get_stylesheet(self):
        """Return filename of CSS stylesheet."""
        dirs = os.getcwd().split('/')
        backsteps=(len(dirs)-dirs.index('astroweb')-1)
        return "../"*backsteps+'css/style.css'

Skeleton.get_stylesheet=get_stylesheet

def get_sitelinks(myfile, rootdir, relthis):
    sitelinks=[]
    # You can add entries via file banner.h
    # if you want
    if os.path.exists('banner.h'):
        banneritems=file('banner.h').readlines()
        for item in banneritems:
            link, name = item.split()
            sitelinks.append((link, name))

    dirs = os.getcwd().split('/')

    whereisit=0
    if 'it' in dirs:
        whereisit=dirs.index('it')
        backstep="../"*(len(dirs)-whereisit)
        new_lang, new_name ='en', '<img height="15" src="%simg/uk.png" alt="English Version"></img>' % backstep
    if 'en' in dirs:
        whereisit=dirs.index('en')
        backstep="../"*(len(dirs)-whereisit)
        new_lang, new_name ='it', '<img height="15" src="%simg/it.png" alt="Versione Italiana"></img>' % backstep
    if whereisit !=0:
        backstep=backstep+new_lang
        for i in range(whereisit+1, len(dirs)): backstep=os.path.join(backstep, dirs[i])
        backstep=os.path.join(backstep, os.path.splitext(myfile)[0]+".html")
        sitelinks.append((backstep, new_name))

    return sitelinks

class Public_ap_generator(Skeleton, Sidebar, Banner):
    AUTHOR = 'Alessandro Pisa'
    EMAIL = 'pisa@fe.infn.it'

    def __init__(self, file, rootdir, relthis):
        root, ext = os.path.splitext(file)
        html = root + '.html'
        p = self.__parser = HTParser(file, self.AUTHOR, self.EMAIL)
        f = self.__linkfixer = LinkFixer(html, rootdir, relthis)
        self.__body = None
        self.__cont = None
        # Calculate the sidebar links, adding a few of our own.
        self.__d = {'rootdir': rootdir}
        dir = os.path.split(file)[0]
        p.process_sidebar()
        p.sidebar.append(BLANKCELL)
        # It is important not to have newlines between the img tag and the end
        # end center tags, otherwise layout gets messed up.
        self.__linkfixer.massage(p.sidebar, self.__d)
        Sidebar.__init__(self, p.sidebar)
        p.sidebar.append(BLANKCELL)
        # Fix up our site links, no relthis because the site links are
        # relative to the root of our web pages.
        sitelinks=get_sitelinks(file, rootdir, relthis)
        sitelink_fixer = LinkFixer(f.myurl(), rootdir)
        sitelink_fixer.massage(sitelinks, self.__d, aboves=1)
        Banner.__init__(self, sitelinks)
        # kludge!
        for i in range(len(p.sidebar)-1, -1, -1):
            if p.sidebar[i] == 'Email Us':
                p.sidebar[i] = 'Email me'
                break
        p.sidebar.append((None, ''))

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
        corner='corner-en.png'
        dirs = os.getcwd().split('/')
        backsteps=(len(dirs)-dirs.index('astroweb')-1)
        if 'it' in dirs: corner='corner-it.png'
        return '''
<center>
    <a class="corner" href="index.html">
    <img class="corner" alt="Astrofisica a Ferrara" border="0"
         src="%s/img/%s"></a></center>''' % ("../"*backsteps,corner)

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
