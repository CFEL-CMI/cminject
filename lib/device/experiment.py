class experiment:
  """This class carry out the info about the experiment you want to simulate. 
     For example, the devices and the interaction fields you will include.
     The output files will have the name of and date of your experiment""" 

  def __init__(self, name, date, devices, fields=None):
    self.name = name
    self.date = date
    self.devices = devices
    self.field = field
