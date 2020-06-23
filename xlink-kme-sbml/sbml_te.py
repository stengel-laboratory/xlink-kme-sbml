# %%
import pandas as pd
import tellurium as te

# %%
exp = 'c3b'

# %%
rr = te.loadSBMLModel(f"model_{exp}_asa.xml")

# %%
with open(
    f"/home/kai/Coding/xlink_kme/xlink_kme/model_{exp}_antimony_asa.txt", "w"
) as f:
    f.write(rr.getAntimony())
# %%
#rr.draw(layout="dot")

# %%
results = rr.simulate(0, 50000000)
print("convergence", max(abs(rr.dv())))

# %%
df_res = pd.DataFrame(results, columns=results.colnames)
# %%
df_res.tail(1).to_csv(f'{exp}_final_frame_asa.csv')



# %%
