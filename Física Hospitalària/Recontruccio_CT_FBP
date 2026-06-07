#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 18:08:27 2026

@author: davidrosamolina
"""

# CT amb imatge DICOM real — Física Hospitalària
# NECESSITEM pip install pydicom scikit-image numpy matplotlib

import numpy as np
import matplotlib.pyplot as plt
import pydicom
import pydicom.data
from skimage.transform import radon, iradon
from skimage.data import shepp_logan_phantom
from skimage.metrics import structural_similarity as ssim
from skimage.metrics import peak_signal_noise_ratio as psnr

# 1. CARREGA IMATGE DICOM REAL (exemple a pydicom)
# Pydicom ens dona una llesca CT real de tòrax de mostra
import os, glob, pydicom
base = os.path.dirname(pydicom.data.__file__)
cts  = glob.glob(os.path.join(base, '**', 'CT_small.dcm'), recursive=True)
ruta = cts[0]
ds   = pydicom.dcmread(ruta)

# Extreu els píxels i aplica la rescala HU (Hounsfield Units)
imatge_raw = ds.pixel_array.astype(float)

# Aplica RescaleSlope i RescaleIntercept si existeixen
slope     = float(getattr(ds, 'RescaleSlope',     1))
intercept = float(getattr(ds, 'RescaleIntercept', 0))
imatge_hu = imatge_raw * slope + intercept

print(f"Pacient: {getattr(ds, 'PatientName', 'Anònim')}")
print(f"Mida imatge: {imatge_hu.shape}")
print(f"Rang HU: {imatge_hu.min():.0f} a {imatge_hu.max():.0f}")
print(f"Gruix de llesca: {getattr(ds, 'SliceThickness', '?')} mm")

# Normalitza a [0, 1] per la reconstrucció
vmin, vmax = -1000, 2000   # rang clínic típic HU
imatge = np.clip((imatge_hu - vmin) / (vmax - vmin), 0, 1)

# Retalla a quadrat (necessari per radon)
n   = min(imatge.shape)
cx  = imatge.shape[0] // 2
cy  = imatge.shape[1] // 2
imatge = imatge[cx - n//2 : cx + n//2,
                cy - n//2 : cy + n//2]


#phantom = shepp_logan_phantom()
#imatge = (phantom - phantom.min()) / (phantom.max() - phantom.min())
#per al phantom matemàtic posar això i treure tot lo del principi
# 2. SINOGRAMA
theta_180 = np.linspace(0., 180., 180, endpoint=False)
sinogram  = radon(imatge, theta=theta_180) #substituir imatge per phantom si volem l'altre model

# 3. RECONSTRUCCIONS FBP
angles = [10, 30, 60, 180]
recons = {}
for n in angles:
    theta_n   = np.linspace(0., 180., n, endpoint=False)
    sino_n    = radon(imatge, theta=theta_n)
    recon     = iradon(sino_n, theta=theta_n, filter_name='ramp')
    recons[n] = np.clip(recon, 0, 1)

# 4. SIMULACIÓ BAIXA DOSI
np.random.seed(42)
recons_noisy = {}
for n in angles:
    noise           = np.random.normal(0, 0.05, recons[n].shape)
    recons_noisy[n] = np.clip(recons[n] + noise, 0, 1)

# 5. MÈTRIQUES 
print("\n" + "=" * 55)
print(f"{'Angles':>8} | {'SSIM neta':>10} | {'SSIM soroll':>11} | {'PSNR neta':>9} | {'PSNR soroll':>10}")
print("-" * 55)
for n in angles:
    s1 = ssim(imatge, recons[n],       data_range=1.0)
    s2 = ssim(imatge, recons_noisy[n], data_range=1.0)
    p1 = psnr(imatge, recons[n],       data_range=1.0)
    p2 = psnr(imatge, recons_noisy[n], data_range=1.0)
    print(f"{n:>8} | {s1:>10.4f} | {s2:>11.4f} | {p1:>8.2f}dB | {p2:>9.2f}dB")
print("=" * 55)

# 6. FIGURA 1: Imatge real a sinograma a FBP
fig, axes = plt.subplots(3, 1, figsize=(4, 13))
fig.suptitle("Phantom, sinograma i reconstrucció",
             fontsize=12, fontweight='bold')

axes[0].imshow(imatge, cmap='gray', vmin=0, vmax=1)
axes[0].set_title("Shepp Logan Phantom", fontsize=10)
axes[0].axis('off')

axes[1].imshow(sinogram, cmap='gray', aspect='auto',
               extent=[0, 180, sinogram.shape[0], 0])
axes[1].set_title("Sinograma\n(180 angles, 0°–180°)", fontsize=10)
axes[1].set_xlabel("Angle (°)")
axes[1].set_ylabel("Posició detector")

axes[2].imshow(recons[180], cmap='gray', vmin=0, vmax=1)
axes[2].set_title("Reconstrucció FBP\n(180 angles, alta dosi)", fontsize=10)
axes[2].axis('off')

plt.tight_layout()
plt.savefig("figura1_recontruccio.png", dpi=300, bbox_inches='tight')
plt.show()

# 7. FIGURA 2: Efecte dels angles amb la imatge real
fig, axes = plt.subplots(2, 4, figsize=(14, 7))
fig.suptitle("Efecte del nombre d'angles a la reconstrucció",
             fontsize=12, fontweight='bold')

for i, n in enumerate(angles):
    s = ssim(imatge, recons[n],       data_range=1.0)
    p = psnr(imatge, recons[n],       data_range=1.0)
    axes[0, i].imshow(recons[n], cmap='gray', vmin=0, vmax=1)
    axes[0, i].set_title(f"{n} angles\nSSIM={s:.3f} PSNR={p:.1f}dB", fontsize=9)
    axes[0, i].axis('off')

    s2 = ssim(imatge, recons_noisy[n], data_range=1.0)
    p2 = psnr(imatge, recons_noisy[n], data_range=1.0)
    axes[1, i].imshow(recons_noisy[n], cmap='gray', vmin=0, vmax=1)
    axes[1, i].set_title(f"+ soroll (baixa dosi)\nSSIM={s2:.3f} PSNR={p2:.1f}dB",
                         fontsize=9)
    axes[1, i].axis('off')

axes[0, 0].set_ylabel("Alta dosi",  fontsize=10, labelpad=10)
axes[1, 0].set_ylabel("Baixa dosi", fontsize=10, labelpad=10)

plt.tight_layout()
plt.savefig("figura2_angles.png", dpi=300, bbox_inches='tight')
plt.show()

# 8. FIGURA 3: Mètriques comparatives
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
fig.suptitle("Mètriques de qualitat imatge CT real",
             fontsize=12, fontweight='bold')

ssim_net   = [ssim(imatge, recons[n],       data_range=1.0) for n in angles]
ssim_noise = [ssim(imatge, recons_noisy[n], data_range=1.0) for n in angles]
psnr_net   = [psnr(imatge, recons[n],       data_range=1.0) for n in angles]
psnr_noise = [psnr(imatge, recons_noisy[n], data_range=1.0) for n in angles]

x, w = np.arange(len(angles)), 0.35

ax1.bar(x - w/2, ssim_net,   w, label='Alta dosi',  color='steelblue')
ax1.bar(x + w/2, ssim_noise, w, label='Baixa dosi', color='salmon')
ax1.set_xticks(x)
ax1.set_xticklabels([f"{n} angles" for n in angles])
ax1.set_ylabel("SSIM (0–1)")
ax1.set_title("SSIM: com més alt, millor")
ax1.set_ylim(0, 1)
ax1.axhline(0.5, color='gray', linestyle='--', linewidth=0.8,
            label='llindar acceptable')
ax1.legend()

ax2.bar(x - w/2, psnr_net,   w, label='Alta dosi',  color='steelblue')
ax2.bar(x + w/2, psnr_noise, w, label='Baixa dosi', color='salmon')
ax2.set_xticks(x)
ax2.set_xticklabels([f"{n} angles" for n in angles])
ax2.set_ylabel("PSNR (dB)")
ax2.set_title("PSNR: com més alt, millor")
ax2.axhline(20, color='gray', linestyle='--', linewidth=0.8,
            label='llindar clínic ~20dB')
ax2.legend()

plt.tight_layout()
plt.savefig("figura3_metriques.png", dpi=300, bbox_inches='tight')
plt.show()

print("\nFigures guardades. Imatge DICOM real de pydicom (CT_small.dcm).")
