# Extract daily data from snapshots taken every 3 hours
for i in range_inc(1, 16):
    with open(f'../data/jpl/testing/hourly/observer-asteroid-{i:03d}-palomar.txt') as fh:
        lines = fh.readlines()
        lines_daily = lines[0:1] + lines[1::8]
        with open(f'../data/jpl/testing/daily/observer-asteroid-{i:03d}-palomar.txt', mode='w') as fh_out:
            fh_out.writelines(lines_daily)