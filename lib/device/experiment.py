class Experiment:
  """This class carry out the info about the experiment you want to simulate. 
     For example, the devices and the interaction fields you will include.
     The output files will have the name of and date of your experiment""" 

  def __init__(self, name, date, detector, source, adl=None, field=None):
    self.name = name
    self.date = date
    self.field = []
    self.device = []

  def add_device(self, device):
    self.device.append(device)

  def add_field(self, device):
    self.field.append(field)

  def run(self):
   pass
