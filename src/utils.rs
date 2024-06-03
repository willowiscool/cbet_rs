/// This function is an implementation of interp that takes a closure which, given an index,
/// returns an element at that index. See the doc for interp for more info.
pub fn interp_closure<F1, F2>(y: F1, ysize: usize, x: F2, xsize: usize, xp: f64) -> f64
where
    F1: Fn(usize) -> f64,
    F2: Fn(usize) -> f64
{
    assert_eq!(xsize, ysize);

    if x(0) <= x(xsize-1) {
        // x monotonically increase
        if xp <= x(0) {
            return y(0);
        } else if xp >= x(xsize-1) {
            return y(ysize-1);
        }

        let mut low = 0;
        let mut high = xsize-1;
        let mut mid = (low+high)>>1;
        while low < high-1 {
            if x(mid) >= xp {
                high = mid;
            } else {
                low = mid;
            }
            mid = (low+high)>>1;
        }
        assert!((xp >= x(mid)) && (xp <= x(mid+1)));
        return y(mid) + ((y(mid+1)-y(mid)) / (x(mid+1)-x(mid)) * (xp-x(mid)));
    } else {
        if xp >= x(0) {
            return y(0);
        } else if xp <= x(xsize-1) {
            return y(ysize-1);
        }

        let mut low = 0;
        let mut high = xsize-1;
        let mut mid = (low+high)>>1;
        while low < high-1 {
            if x(mid) <= xp {
                low = mid;
            } else {
                high = mid;
            }
            mid = (low + high) >> 1;
        }
        assert!((xp <= x(mid)) && (xp >= x(mid+1)));
        return y(mid) + ((y(mid+1)-y(mid))/(x(mid+1)-x(mid))*(xp-x(mid)));
    }
}

/// This function is copied from c++ impl., may be able to be made nicer or eliminated entirely.
///
/// Assumes x either monotonically increases or decreases. y, x is a set of points. Returns the
/// interpolated y value for the x value given.
///
/// NOTE: does not handle errors. panics instead. (TODO fix??)
/// 
/// Possible optimization (as mentioned in ray_trace.rs, before interp is first used): let
/// interp take closures for y, x so that if, for example, your y/x are a column of a vector,
/// you don't have to copy ALL the values beforehand.
pub fn interp(y: &Vec<f64>, x: &Vec<f64>, xp: f64) -> f64 {
    assert_eq!(x.len(), y.len());

    if x[0] <= *x.last().unwrap() {
        // x monotonically increase
        if xp <= x[0] {
            return y[0];
        } else if xp >= *x.last().unwrap() {
            return *y.last().unwrap();
        }

        let mut low = 0;
        let mut high = x.len() - 1;
        let mut mid = (low + high) >> 1;
        while low < high - 1 {
            if x[mid] >= xp {
                high = mid;
            } else {
                low = mid;
            }
            mid = (low + high) >> 1;
        }

        assert!((xp >= x[mid]) && (xp <= x[mid + 1]));
        return y[mid] + ((y[mid+1] - y[mid]) / (x[mid+1] - x[mid]) * (xp - x[mid]));
    } else {
        if xp >= x[0] {
            return y[0];
        } else if xp <= *x.last().unwrap() {
            return *y.last().unwrap();
        }

        let mut low = 0;
        let mut high = x.len() - 1;
        let mut mid = (low + high) >> 1;
        while low < high - 1 {
            if x[mid] <= xp {
                low = mid;
            } else {
                high = mid;
            }
            mid = (low + high) >> 1;
        }

        assert!((xp <= x[mid]) && (xp >= x[mid+1]));
        return y[mid] + ((y[mid+1]-y[mid])/(x[mid+1]-x[mid])*(xp-x[mid]));
    }
}
