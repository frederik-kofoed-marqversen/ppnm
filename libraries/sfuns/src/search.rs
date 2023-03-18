pub fn binary_search<T: PartialOrd>(sorted_array: &Vec<T>, item: T) -> Result<usize, ()> {
    let (mut left, mut right) = (0, sorted_array.len() - 1);
    while left != right {
        let midpoint: usize = (left + right + 1) / 2;  // this is ceil((left+right)/2)
        if sorted_array[midpoint] > item {
            right = midpoint - 1;
        } else {
            left = midpoint;
        }
    }
    if sorted_array[left] == item {return Ok(left);}
    return Err(());
}

pub fn binary_search_leftmost<T: PartialOrd>(sorted_array: &Vec<T>, item: T) -> usize {
    // returns the number of elements to the left of (less than) 'item'
    let (mut left, mut right) = (0, sorted_array.len());
    
    while left < right {
        let midpoint: usize = (right + left) / 2;
        if sorted_array[midpoint] < item {
            left = midpoint + 1;
        } else {
            right = midpoint;
	    }
    }
    return left
}

pub fn binary_search_rightmost<T: PartialOrd>(sorted_array: &Vec<T>, item: T) -> usize {
    let (mut left, mut right) = (0, sorted_array.len());
    
    while left < right {
        let midpoint: usize = (right + left) / 2;
        if sorted_array[midpoint] < item {
            right = midpoint;
        } else {
            left = midpoint + 1;
	    }
    }
    return right - 1
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_binary_search_1() {
        let array: Vec<u64> = (4..20).collect();
        assert_eq!(
            binary_search(&array, 18),
            Ok(14)
        );
    }

    #[test]
    fn test_binary_search_2() {
        let array: Vec<u64> = (4..20).collect();
        assert_eq!(
            binary_search(&array, 22),
            Err(())
        );
    }

    #[test]
    fn test_binary_search_3() {
        let array: Vec<u64> = (4..20).collect();
        assert_eq!(
            binary_search(&array, 0),
            Err(())
        );
    }

    #[test]
    fn test_binary_search_leftmost_1() {
        let array = vec![2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        assert_eq!(
            binary_search_leftmost(&array, 3.0),
            1
        );
    }

    #[test]
    fn test_binary_search_leftmost_2() {
        let array = vec![2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        assert_eq!(
            binary_search_leftmost(&array, 3.2),
            2
        );
    }

    #[test]
    fn test_binary_search_leftmost_3() {
        let array = vec![2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        assert_eq!(
            binary_search_leftmost(&array, 0.0),
            0
        );
    }

    #[test]
    fn test_binary_search_leftmost_4() {
        let array = vec![2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        assert_eq!(
            binary_search_leftmost(&array, 2.0),
            0
        );
    }
}