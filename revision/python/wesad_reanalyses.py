# ============================================================
# WESAD reanalyses from SAVED canonical Model-D fit outputs
# (post-processing only: muhat/sigmahat are the fitted values;
#  no refitting is done here). Produces real numbers for the
# reviewer response: endpoint counts, participant-clustered
# bootstrap AUC CIs, threshold-leakage sensitivity, Model-D
# calibration for the exceedance event, 77-feature list.
# ============================================================
import pandas as pd, numpy as np, json, os
from scipy.stats import norm
from sklearn.metrics import roc_auc_score, brier_score_loss

rng = np.random.default_rng(2026)
W = "revision/data"
OUT = "revision/outputs"
os.makedirs(OUT, exist_ok=True)

rs  = pd.read_csv(f"{W}/risk_scores.csv")
lab = pd.read_csv(f"{W}/label_tbl.csv")
datfit = pd.read_csv(f"{W}/dat_fit.csv")

# ---- merge saved fitted mean/sd with labels ----
m = rs.merge(lab, on=["s","t"], how="inner").copy()
m["grp"] = (m["label"]==2).astype(int)          # stress = positive
# subject-specific 70th-percentile threshold over FULL timeframe (as in paper)
m["cs70_full"] = m.groupby("s")["y"].transform(lambda v: v.quantile(0.70))
m["risk70"] = 1 - norm.cdf((m["cs70_full"] - m["muhat"]) / m["sigmahat"])

def pooled_auc(y, s):   # fixed orientation: higher score -> stress
    return roc_auc_score(y, s)

# orientation check
a_risk = pooled_auc(m["grp"], m["risk70"]); a_hr = pooled_auc(m["grp"], m["y"])
print(f"[orientation] AUC risk70={a_risk:.3f}  AUC HR={a_hr:.3f}  (both should be >0.5)")

# =====================================================================
# (1) ENDPOINT DEFINITION + COUNTS BY PARTICIPANT  (Reviewer 3.1)
# =====================================================================
lab_names = {0:"transition/other",1:"baseline",2:"stress",3:"amusement"}
cnt = (m.groupby(["s","label"]).size().unstack(fill_value=0)
         .rename(columns=lab_names))
for c in ["baseline","stress","amusement","transition/other"]:
    if c not in cnt.columns: cnt[c]=0
cnt = cnt[["baseline","stress","amusement","transition/other"]]
cnt["total_windows"] = cnt.sum(axis=1)
cnt["stress(pos)"]   = cnt["stress"]
cnt["non-stress(neg)"] = cnt["total_windows"] - cnt["stress"]
cnt.loc["TOTAL"] = cnt.sum(axis=0)
cnt.to_csv(f"{OUT}/counts_by_participant.csv")
print("\n=== (1) Window/class counts by participant ===")
print(cnt.to_string())

# =====================================================================
# (2) PARTICIPANT-CLUSTERED BOOTSTRAP AUC CIs  (Reviewer 3.1)
# =====================================================================
subs = sorted(m["s"].unique())
def cluster_boot(score_col, B=3000):
    aucs=[]
    grp_by = {s: m.loc[m.s==s,"grp"].values for s in subs}
    scr_by = {s: m.loc[m.s==s,score_col].values for s in subs}
    for _ in range(B):
        pick = rng.choice(subs, size=len(subs), replace=True)
        y = np.concatenate([grp_by[s] for s in pick])
        sc= np.concatenate([scr_by[s] for s in pick])
        if y.sum()==0 or y.sum()==len(y): continue
        aucs.append(roc_auc_score(y,sc))
    return np.percentile(aucs,[2.5,50,97.5]), np.std(aucs)

def cluster_boot_delta(B=3000):
    d=[]
    g={s:m.loc[m.s==s,"grp"].values for s in subs}
    r={s:m.loc[m.s==s,"risk70"].values for s in subs}
    h={s:m.loc[m.s==s,"y"].values for s in subs}
    for _ in range(B):
        pick=rng.choice(subs,size=len(subs),replace=True)
        y=np.concatenate([g[s] for s in pick])
        if y.sum()==0 or y.sum()==len(y): continue
        rr=np.concatenate([r[s] for s in pick]); hh=np.concatenate([h[s] for s in pick])
        d.append(roc_auc_score(y,rr)-roc_auc_score(y,hh))
    return np.percentile(d,[2.5,50,97.5])

ci_risk,_ = cluster_boot("risk70")
ci_hr,_   = cluster_boot("y")
ci_d      = cluster_boot_delta()
print("\n=== (2) Participant-clustered bootstrap 95% CIs (3000 reps, resampling 15 subjects) ===")
print(f"AUC risk70 : point={a_risk:.3f}  95%CI=[{ci_risk[0]:.3f}, {ci_risk[2]:.3f}]")
print(f"AUC HR     : point={a_hr:.3f}  95%CI=[{ci_hr[0]:.3f}, {ci_hr[2]:.3f}]")
print(f"Delta(risk-HR): point={a_risk-a_hr:.3f}  95%CI=[{ci_d[0]:.3f}, {ci_d[2]:.3f}]")
boot = {"auc_risk":a_risk,"ci_risk":list(ci_risk),"auc_hr":a_hr,"ci_hr":list(ci_hr),
        "delta":a_risk-a_hr,"ci_delta":list(ci_d)}
json.dump(boot, open(f"{OUT}/bootstrap_ci.json","w"), indent=2)

# =====================================================================
# (3) THRESHOLD-LEAKAGE SENSITIVITY  (Reviewer 2.2)
#     Recompute c_s three ways; muhat/sigmahat unchanged.
# =====================================================================
print("\n=== (3) Threshold choice sensitivity (information-leak check) ===")
res=[]
# (a) full-timeframe quantile (paper)
for p in [0.70]:
    m["cs"]=m.groupby("s")["y"].transform(lambda v: v.quantile(p))
    m["rk"]=1-norm.cdf((m["cs"]-m["muhat"])/m["sigmahat"])
    res.append(("full-timeframe q0.70", roc_auc_score(m["grp"],m["rk"])))
# (b) burn-in only: first K windows per subject (protocol starts at baseline)
for frac in [0.20, 0.30]:
    def thr_burn(v, frac=frac):
        k=max(5,int(np.floor(frac*len(v))))
        return np.quantile(v.iloc[:k], 0.70)
    cs = m.groupby("s")["y"].transform(lambda v: thr_burn(v))
    rk = 1-norm.cdf((cs-m["muhat"])/m["sigmahat"])
    res.append((f"burn-in first {int(frac*100)}% q0.70", roc_auc_score(m["grp"],rk)))
# (c) baseline-window quantile (label==1 windows only) as a clean non-stress reference
def thr_base(sdf):
    base = sdf.loc[sdf.label==1,"y"]
    ref  = base if len(base)>=5 else sdf["y"]
    return np.quantile(ref,0.70)
csb = m.groupby("s",group_keys=False).apply(lambda d: pd.Series(thr_base(d), index=d.index))
rkb = 1-norm.cdf((csb-m["muhat"])/m["sigmahat"])
res.append(("baseline-window q0.70", roc_auc_score(m["grp"],rkb)))
for name,auc in res: print(f"  c_s from {name:32s}: AUC(stress) = {auc:.3f}")
pd.DataFrame(res,columns=["threshold_rule","AUC_stress"]).to_csv(f"{OUT}/threshold_sensitivity.csv",index=False)

# =====================================================================
# (4) MODEL-D CALIBRATION for the actual exceedance event  (Reviewer 3.3 template)
#     event = 1{ y > c_s(0.70) };  risk70 predicts it.
# =====================================================================
m["event70"] = (m["y"] > m["cs70_full"]).astype(int)
brier = brier_score_loss(m["event70"], m["risk70"])
# reliability by decile
m["bin"]=pd.qcut(m["risk70"],10,labels=False,duplicates="drop")
rel=m.groupby("bin").agg(mean_pred=("risk70","mean"),obs_rate=("event70","mean"),n=("event70","size"))
auc_event = roc_auc_score(m["event70"], m["risk70"])
print("\n=== (4) Model D calibration for exceedance event y>c_s(0.70) ===")
print(f"  Brier score = {brier:.4f} | AUC(event) = {auc_event:.3f} | event rate = {m['event70'].mean():.3f}")
rel.to_csv(f"{OUT}/modelD_calibration_reliability.csv")
json.dump({"brier":brier,"auc_event":auc_event,"event_rate":float(m["event70"].mean())},
          open(f"{OUT}/modelD_calibration.json","w"),indent=2)

# =====================================================================
# (5) 77-FEATURE LIST  (Reviewer 3.4)
# =====================================================================
feats=[c for c in datfit.columns if c not in ("s","t","y")]
def block(f):
    u=f.upper()
    if u.startswith("HRV") or u.startswith("ECG"): return "ECG / HRV"
    if u.startswith(("RSP","BRV","INSP","EXP","DUTY")): return "Respiration & BRV"
    if u.startswith(("EDA","SCL","SCR")): return "EDA"
    if u.startswith("EMG"): return "EMG"
    if u.startswith("ACC"): return "Accelerometry"
    if u.startswith("TEMP") or "TEMP" in u and u.startswith("TEMP"): return "Skin temperature"
    if u.startswith("CORR"): return "Cross-channel correlations"
    if "TEMP" in u: return "Skin temperature"
    return "Other"
fl=pd.DataFrame({"feature":feats}); fl["block"]=fl["feature"].map(block)
fl=fl.sort_values(["block","feature"])
fl.to_csv(f"{OUT}/feature_list_77.csv",index=False)
print("\n=== (5) 77 features by block ===")
print(fl["block"].value_counts().to_string())
print("total features:",len(feats))
print("\nDONE. Outputs in", OUT)
