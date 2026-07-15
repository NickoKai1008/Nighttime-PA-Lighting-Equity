from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent

SCRIPTS = (
    Path("Fig1_gender_PA_disparities/scripts/gender_footprint.py"),
    Path("Fig2_nighttime_gender_gap/scripts/plot_gender_inequality_quadrant.py"),
    Path("Fig2_nighttime_gender_gap/scripts/plot_lorenz_curves.py"),
    Path("Fig3_nonlinear_response/scripts/fig3_nonlinear_response.py"),
    Path("Fig4_overlay_forest/scripts/overlay_forest.py"),
    Path("Fig5_ALAN_reallocation_curve_actionmap/scripts/fig5_ALAN_reallocation.py"),
)


def main() -> None:
    for relative_path in SCRIPTS:
        script = ROOT / relative_path
        print(f"Running {relative_path.as_posix()}", flush=True)
        subprocess.run(
            [sys.executable, str(script)],
            cwd=script.parents[1],
            check=True,
        )

    print("All figure scripts completed successfully.")


if __name__ == "__main__":
    main()
