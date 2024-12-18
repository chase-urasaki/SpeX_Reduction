# #%%
# import numpy as np 

# #order1 = np.asarray([x1:=380, x2:=500, y1:=460, y2:=1000]) 
# class order_traces:

# order1 = np.asarray([380, 500, 460, 1000]) 
# # %%
# # %%
#%%
import numpy as np

class OrderTraces:
    def __init__(self):
        # Define order1 as a NumPy array
        self.order1 = np.asarray([x1 := 380, x2 := 500, y1 := 460, y2 := 1000])
    
    def get_order(self, name):
        """Retrieve an order trace by name."""
        if hasattr(self, name):
            return getattr(self, name)
        else:
            raise AttributeError(f"No such order trace: {name}")

# %%
