Mon 13 Jun 2016 04:41:17 PM EDT

* Also, it makes sense to do the camera pointings by __orbit__, not by "sector". The relevant
factor of 2 is line 116 and 119 of `eclip_observe`.
  This changes the NUMBER by a factor of two, and also means we'll rewrite all the input
lists to be by orbit instead of by sector. The reason to do this is to have something consistent
with \hemis14d.
  As far as how this affects our shoddy job at dropping fields from the earth/moon crossings,
let's just say it won't. We already decided that it would be the worst-affected SECTORS, so we
can just cascade that down to the orbit-specific pointings on a per-camera basis.

If we wanted to be more precise, we would probably only drop like 1 in 2 of the orbits
per sector, since that seems on average closer to the fraction of dropped time than 0 or 100%.

OFC, our data from Jacobi has this "300ct/s/px" cutoff, which underestimates the fraction of
time that these crossings will be causing bad effects on the photometry.

Here's what we do for the memo:
  * Make pointings `*_coord.dat files` for all 6 selected 1yr scenarios, without earth/moon drops.
    Call them `*_orbits.dat` files.
    X shemi_nhemi_orbits.dat
    X nhemi_orbits.dat
    X hemis14d_orbits.dat
    X npole_orbits.dat
    X shemiAvoid
    X elong: note there's ambiguity here inre: when the "good 7 sectors" are during the year.
      Since I don't have clarity from jacobi on the specific orbital phasing, and don't particularly
      care (it only affects continuity in the 6 npole sectors, which are less the point)
      So let's just do 6 npole sectors first, then 7 elong.
    X eshort
  * Drop earth moon sectors according to following procedure:
    Anywhere I say `table', I'm referring to `outage_time_table.xlsx`.
    When I'm referring to the outage plots, those are in `EarthMoonSummary300.xlsx`.
    The relevant plots I'm looking at when determing "which sectors have the worst earth/moon
    coverage" are in 
      * Take cumulative number of "lost sectors" from table, look at plot from table to see where
        they are grouped in time (n.b. there's some phase difference btwn the planet detection sim
        and the orbit sim, but it doesn't matter since it's fixed)
      * For each camera, drop the fields that have the worst crossings (i.e. those that correspond
        to the largest times of fractional outage, per sector). 
        The total number of fields dropped should correspond to the cumulative number above 
        (i.e. if "2 sectors" cumulative, drop 4 orbit camera pointings).
    These will be called `*_orbout.dat`, since these are orbits with outages accounted for.
    X shemi_nhemi
    X nhemi
    X npole
    X shemiAvoid: chose to drop the sequential fields early on (& not opposite ends of year)
    X elong: recall pole continuous for first 6 sectors. To avoid latitude bias, just drop 4
      of the camera 0 fields across every-other sector.
    X eshort: same as elong, but only 1 camera 1 sector, and 1 camera 2 sector.
    X hemis14d: drop the expected fields early in the mission.

