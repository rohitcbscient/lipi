import pyvo as vo
import pandas as pd
from datetime import datetime, timedelta

# --------------------------------------------------
# Connect to the MWA TAP service
# --------------------------------------------------
mwa_tap_service = vo.dal.TAPService(
    "https://vo.mwatelescope.org/mwa_asvo/tap"
)

# --------------------------------------------------
# Time windows (UT)
# --------------------------------------------------
windows = [
    ("01:00:00", "02:00:00"),
    ("03:00:00", "04:00:00"),
    ("06:00:00", "07:00:00"),
]

# --------------------------------------------------
# Loop through all days in 2014
# --------------------------------------------------
start_date = datetime(2014, 1, 1)
end_date = datetime(2015, 1, 1)

all_rows = []

current_date = start_date

while current_date < end_date:

    date_str = current_date.strftime("%Y-%m-%d")

    for start_ut, stop_ut in windows:

        t0 = f"{date_str} {start_ut}"
        t1 = f"{date_str} {stop_ut}"

        query = f"""
        SELECT
            obs_id,
            obsname,
            starttime_utc,
            stoptime_utc,
            projectid,
            calibration,
            calibrators,
            deleted_flag
        FROM mwa.observation
        WHERE starttime_utc >= '{t0}'
        AND starttime_utc <  '{t1}'
        AND projectid = 'G0002'
        AND deleted_flag = '0'
        ORDER BY starttime_utc
        """

        try:
            results = mwa_tap_service.search(query)

            # Keep only non-empty results
            if len(results) > 0:
                for row in results:
                    all_rows.append(dict(row))

                print(
                    f"{date_str} {start_ut}-{stop_ut}: "
                    f"{len(results)} observations"
                )

        except Exception as e:
            print(
                f"Failed for {date_str} "
                f"{start_ut}-{stop_ut}: {e}"
            )

    current_date += timedelta(days=1)

# --------------------------------------------------
# Save results
# --------------------------------------------------
if len(all_rows) > 0:

    df = pd.DataFrame(all_rows)

    # Remove duplicate observations if any
    df = df.drop_duplicates(subset=["obs_id"])

    # Sort chronologically
    df = df.sort_values("starttime_utc")

    output_file = "MWA_G0002_2014_solar_obs_metadata.csv"
    df.to_csv(output_file, index=False)

    print(f"\nSaved {len(df)} observations to:")
    print(output_file)

else:
    print("No matching observations found.")