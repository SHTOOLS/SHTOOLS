"""
Date and time utils for SHTOOLS
"""
import datetime as _datetime


def _yyyymmdd_to_year_fraction(date):
    """Convert YYYMMDD.DD date string or float to YYYY.YYY"""
    date = str(date)
    if '.' in date:
        date, residual = str(date).split('.')
        residual = float('0.' + residual)
    else:
        residual = 0.0

    date = _datetime.datetime.strptime(date, '%Y%m%d')
    date += _datetime.timedelta(days=residual)

    year = date.year
    year_start = _datetime.datetime(year=year, month=1, day=1)
    next_year_start = _datetime.datetime(year=year + 1, month=1, day=1)
    year_duration = next_year_start - year_start

    year_elapsed = date - year_start
    fraction = year_elapsed / year_duration

    return year + fraction
