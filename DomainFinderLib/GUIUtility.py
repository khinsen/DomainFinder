import Tkinter, Dialog, tkFileDialog
import os, string

GUIError = 'GUIError'


class FilenameEntry(Tkinter.Frame):

    def __init__(self, master, text, browse_pattern = '*', must_exist = 1):
        self.pattern = browse_pattern
        self.must_exist = must_exist
        Tkinter.Frame.__init__(self, master)
        Tkinter.Label(self, text=text).pack(side=Tkinter.LEFT)
        self.filename = Tkinter.StringVar()
        Tkinter.Button(self, text="Browse...",
                       command=self.browse).pack(side=Tkinter.RIGHT)
        entry = Tkinter.Entry(self, textvariable=self.filename)
        entry.pack(side=Tkinter.RIGHT, expand=1, fill=Tkinter.X)
        entry.icursor("end")

    def browse(self):
        open = tkFileDialog.Open(self,
                                 filetypes=[("DomainFinder files", self.pattern),
                                             ("All files", "*")],
                                 title = "Choose output file")
        file = open.show()
        if file:
            self.filename.set(file)

    def get(self):
        filename =  self.filename.get()
        if self.must_exist and not os.path.exists(filename):
            Dialog.Dialog(self, title='File not found',
                          text='The file "' + filename + '" does not exist.',
                          bitmap='warning', default=0,
                          strings = ('Cancel',))
            raise ValueError
        return filename


class FloatEntry(Tkinter.Frame):

    def __init__(self, master, text, init = None, lower=None, upper=None,
                 name = None):
        self.text = text
        self.lower = lower
        self.upper = upper
        if name is None:
            name = text
        self.name = name
        if len(text) > 0:
            label_text = string.upper(text[0])+text[1:]+':'
        else:
            label_text = ''
        Tkinter.Frame.__init__(self, master)
        Tkinter.Label(self, text=label_text).pack(side=Tkinter.LEFT)
        self.value = Tkinter.DoubleVar()
        if init is not None:
            self.value.set(init)
        self.entry = Tkinter.Entry(self, textvariable=self.value)
        self.entry.pack(side=Tkinter.RIGHT, anchor=Tkinter.E,
                        expand=1, fill=Tkinter.X)
        self.entry.icursor("end")

    def bind(self, sequence=None, func=None, add=None):
        self.entry.bind(sequence, func, add)

    def set(self, value):
        return self.value.set(value)

    def get(self):
        try:
            value = self.value.get()
        except Tkinter.TclError:
            Dialog.Dialog(self, title='Illegal value',
                          text='The value of "' + self.name +
                               '" must be a number.',
                          bitmap='warning', default=0,
                          strings = ('Cancel',))
            raise ValueError
        range_check = 0
        if self.lower is not None and value < self.lower:
            range_check = -1
        if self.upper is not None and value > self.upper:
            range_check = 1
        if range_check != 0:
            text = 'The value of "' + self.name + '" must not be '
            if range_check < 0:
                text = text + 'smaller than ' + `self.lower` + '.'
            else:
                text = text + 'larger than ' + `self.upper` + '.'
            Dialog.Dialog(self, title='Value out of range', text=text,
                          bitmap='warning', default=0,
                          strings = ('Cancel',))
            raise ValueError
        return value


class IntEntry(FloatEntry):

    def get(self):
        value = FloatEntry.get(self)
        ivalue = int(value)
        if ivalue != value:
            Dialog.Dialog(self, title='Illegal value',
                          text='The value of "' + self.name +
                               '" must be an integer.',
                          bitmap='warning', default=0,
                          strings = ('Cancel',))
            raise ValueError
        return ivalue

class ButtonBar(Tkinter.Frame):

    def __init__(self, master, left_button_list, right_button_list):
        Tkinter.Frame.__init__(self, master, bd=2, relief=Tkinter.SUNKEN)
        for button, action in left_button_list:
            Tkinter.Button(self, text=button,
                           command=action).pack(side=Tkinter.LEFT)
        for button, action in right_button_list:
            Tkinter.Button(self, text=button,
                           command=action).pack(side=Tkinter.RIGHT)


class StatusBar(Tkinter.Frame):

    def __init__(self, master):
        Tkinter.Frame.__init__(self, master, bd=2, relief=Tkinter.RAISED)
        self.text = Tkinter.Label(self, text='')
        self.text.pack(side=Tkinter.LEFT, expand=Tkinter.YES)

    def set(self, text):
        self.text.configure(text = text)
        self.text.update_idletasks()
        self.master.configure(cursor='watch')
        self.update()
        self.update_idletasks()

    def clear(self):
        self.text.configure(text = '')
        self.text.update_idletasks()
        self.master.configure(cursor='top_left_arrow')
        self.update_idletasks()
