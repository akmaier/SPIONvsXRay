"""SPIONvsXRay — simulating CT visibility of iron-oxide nanoparticles.

See SPEC.md for the full experimental plan. Modules:
  config       — all experiment parameters (single source of truth)
  materials    — X-ray materials (soft tissue, bone, iron-loaded tumor)
  phantom      — rabbit-scale digital phantom with tumor + bone
  simulate     — polychromatic forward projection (500-view C-arm)
  reconstruct  — cone-beam FDK
  analyze      — ROI HU/CNR, detectability
  build_dashboard — regenerate docs/ from results/
"""
__version__ = "0.1.0"
