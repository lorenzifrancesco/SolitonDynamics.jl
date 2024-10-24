from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import time
import matplotlib.pyplot as plt
import numpy as np

def exec(path):
  print("hit")
  plt.figure(figsize=(3.5, 2.7))
  plt.savefig("media/shape/soliton.pdf")



###################################
###################################

class FileChangeHandler(FileSystemEventHandler):
    def __init__(self, callback):
        self.callback = callback

    def on_modified(self, event):
        if event.is_directory:
            return
           
        print("="*50)
        print(f"File {event.src_path} has been modified.")
        self.callback(event.src_path)

if __name__ == "__main__":
    exec("")
    event_handler = FileChangeHandler(exec)
    observer = Observer()
    observer.schedule(event_handler, path='/home/lorenzi/sw/SolitonDynamics.jl/results', recursive=False)
    observer.start()

    try:
        while True:
            time.sleep(0.2)
    except KeyboardInterrupt:   
        observer.stop()

    observer.join()