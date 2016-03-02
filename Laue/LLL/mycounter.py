"""Progressbar class provided with a kdialog progressbar
"""
from threading import Thread
from os import popen
from time import sleep

class Progressbar(Thread):
    '''Classed based on Thread from threading

    index: the initial counter value
    update_time: self explaining
    steps: number of steps to do
    title: window title
    name: name given to the progress bar
    '''
    def __init__(self,
                 index=0,
                 update_time=1., # seconds
                 title='mytitle',
                 name='counter',
                 steps=1.e6,
                 CancelButton=True,
                ):
        Thread.__init__(self)
        self.index=index
        self.name=name
        self.title=title
        self.update_time=update_time
        self.steps=steps
        self.stop=False
        self.dcop_id=self.dcop_start()
        if CancelButton: self.showCancelButton()
        self.setTotalSteps()
        self.start()

    def dcop_start(self):
        '''Gets the id from the progress dialog.
        The strings returned by the pipe is something like
        "DCOPRef(kdialog-NNNNN,ProgressDialog)\n"
        '''
        return popen('/usr/bin/kdialog --nograb --geometry 256 --progressbar "%s" --title "%s"' % (self.name, self.title)).read()[:-1]

    def setTotalSteps(self):
        '''Sets the maximum number of steps in the dialog.'''
        popen('/usr/bin/dcop "%s" setTotalSteps %s' % (self.dcop_id, self.steps))

    def setProgress(self):
        '''Updates the dialog with the counted steps'''
        popen('/usr/bin/dcop "%s" setProgress %s' % (self.dcop_id, self.index))

    def showCancelButton(self):
        '''Show the Cancel Button on the dialog'''
        popen('/usr/bin/dcop "%s" showCancelButton true' % self.dcop_id)

    def wasCancelled(self):
        '''Polls dcop to see if the Cancel Button (if it exists...) was pressed'''
        if popen('/usr/bin/dcop "%s" wasCancelled' % self.dcop_id).read() == "true\n": self.stop=True

    def close(self):
        '''Close the dialog'''
        self.stop=True
        popen('/usr/bin/dcop "%s" close' % self.dcop_id)

    def run(self):
        '''Override the run function from Thread class, called by self.start().'''
        while self.stop is False:
            self.setProgress()
            self.wasCancelled()
            sleep(self.update_time)

def test():
    '''Test function used to show the capability of the implemented progress bar object.'''
    pbar=Progressbar(title='test title',
                     name='test window',
                     update_time=.1,
                     steps=1000,
                     )
    for pbar.index in range(1000):
        sleep(.01)
        if pbar.stop: break
    pbar.close()

if __name__=="__main__":
    test()
