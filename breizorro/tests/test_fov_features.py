"""Tests for field-of-view and radial cutoff features in breizorro"""
import numpy as np
import pytest
from breizorro.utils import match_mask_shape, parse_fov, apply_fov_crop, apply_radial_cutoff


class TestMatchMaskShape:
    """Test cases for match_mask_shape function"""

    def test_crop_larger_mask_to_smaller_target(self):
        """Test cropping a larger mask to a smaller target shape"""
        large_mask = np.ones((100, 100))
        target_shape = (50, 50)
        result = match_mask_shape(large_mask, target_shape)
        
        assert result.shape == target_shape, "Shape mismatch!"
        assert np.array_equal(result, np.ones(target_shape)), "Values incorrect!"

    def test_extend_smaller_mask_to_larger_target(self):
        """Test extending a smaller mask to a larger target with zero-padding"""
        small_mask = np.ones((30, 30))
        target_shape = (80, 80)
        result = match_mask_shape(small_mask, target_shape)
        
        assert result.shape == target_shape, "Shape mismatch!"
        assert np.sum(result) == 30*30, "Extended mask should have zeros padding!"
        # Verify the original data is in the top-left corner
        assert np.array_equal(result[:30, :30], small_mask), "Data placement incorrect!"

    def test_same_shape_preservation(self):
        """Test that same-sized masks are preserved"""
        same_mask = np.ones((50, 50))
        target_shape = (50, 50)
        result = match_mask_shape(same_mask, target_shape)
        
        assert result.shape == target_shape
        assert np.array_equal(result, same_mask), "Should be identical!"

    def test_match_with_different_values(self):
        """Test matching with non-uniform mask values"""
        mask = np.arange(100).reshape(10, 10)
        target_shape = (15, 12)
        result = match_mask_shape(mask, target_shape)
        
        assert result.shape == target_shape
        assert np.array_equal(result[:10, :10], mask), "Original data should be preserved!"


class TestParseFOV:
    """Test cases for parse_fov function"""

    def test_valid_fov_string(self):
        """Test parsing a valid FOV string"""
        fov_str = "10,20,90,80"
        result = parse_fov(fov_str)
        
        assert result == (10, 90, 20, 80), "FOV parsing failed!"

    def test_none_fov(self):
        """Test that None FOV returns None"""
        result = parse_fov(None)
        assert result is None, "None should return None!"

    def test_empty_string_fov(self):
        """Test that empty string FOV returns None"""
        result = parse_fov("")
        assert result is None, "Empty string should return None!"

    def test_invalid_fov_wrong_length(self):
        """Test that wrong number of values raises error"""
        with pytest.raises(ValueError, match="exactly 4 comma-separated values"):
            parse_fov("10,20,90")

    def test_invalid_fov_not_string(self):
        """Test that non-string type raises error"""
        with pytest.raises(ValueError, match="must be a string"):
            parse_fov([10, 20, 90, 80])

    def test_invalid_fov_negative_values(self):
        """Test that negative values raise error"""
        with pytest.raises(ValueError, match="Invalid FOV range"):
            parse_fov("-10,20,90,80")

    def test_invalid_fov_reversed_x_range(self):
        """Test that x1 <= x0 raises error"""
        with pytest.raises(ValueError, match="Invalid FOV range"):
            parse_fov("90,20,10,80")

    def test_invalid_fov_reversed_y_range(self):
        """Test that y1 <= y0 raises error"""
        with pytest.raises(ValueError, match="Invalid FOV range"):
            parse_fov("10,80,90,20")


class TestApplyFOVCrop:
    """Test cases for apply_fov_crop function"""

    def test_apply_fov_crop(self):
        """Test basic FOV cropping"""
        test_mask = np.arange(100*100).reshape(100, 100)
        fov_str = "10,20,50,60"
        result, actual_fov = apply_fov_crop(test_mask, fov_str)
        
        assert result.shape == (40, 40), "Cropped shape incorrect!"
        assert np.array_equal(result, test_mask[20:60, 10:50]), "Cropped values incorrect!"
        assert actual_fov == (10, 50, 20, 60), "Actual FOV mismatch!"

    def test_no_fov_applied(self):
        """Test that no FOV returns original mask"""
        test_mask = np.arange(100*100).reshape(100, 100)
        result, actual_fov = apply_fov_crop(test_mask, None)
        
        assert np.array_equal(result, test_mask), "Should return original!"
        assert actual_fov is None, "actual_fov should be None!"

    def test_fov_exceeds_bounds_clamped(self):
        """Test that FOV exceeding bounds is clamped to image dimensions"""
        test_mask = np.ones((100, 100))
        result, actual_fov = apply_fov_crop(test_mask, "10,20,200,150")
        
        assert actual_fov == (10, 100, 20, 100), "Should clamp to max dimensions!"
        assert result.shape == (80, 90), "Clamped crop shape incorrect!"

    def test_fov_completely_outside_image(self):
        """Test that FOV completely outside image raises error"""
        test_mask = np.ones((100, 100))
        
        with pytest.raises(ValueError, match="invalid crop"):
            apply_fov_crop(test_mask, "150,150,200,200")

    def test_fov_partial_outside_clamped(self):
        """Test that partially outside FOV is clamped gracefully"""
        test_mask = np.ones((100, 100))
        result, actual_fov = apply_fov_crop(test_mask, "50,0,150,100")
        
        assert actual_fov == (50, 100, 0, 100), "Should clamp x1 to image width!"
        assert result.shape == (100, 50), "Result shape should match clamped region!"

    def test_fov_full_image(self):
        """Test FOV covering the entire image"""
        test_mask = np.ones((100, 100))
        result, actual_fov = apply_fov_crop(test_mask, "0,0,100,100")
        
        assert np.array_equal(result, test_mask), "Should return full image!"
        assert actual_fov == (0, 100, 0, 100), "Actual FOV should be full image!"


class TestApplyRadialCutoff:
    """Test cases for apply_radial_cutoff function"""

    def test_radial_cutoff_basic(self):
        """Test basic radial cutoff application"""
        test_mask = np.ones((100, 100))
        radius = 30
        result = apply_radial_cutoff(test_mask, radius)
        
        assert result.shape == test_mask.shape, "Shape should be preserved!"
        
        # Check that center is still 1 and corners are 0
        center_y, center_x = 50, 50
        assert result[center_y, center_x] == 1.0, "Center should be 1!"
        assert result[0, 0] == 0.0, "Corner should be 0!"
        
        # Count non-zero pixels (should be approximately pi*r^2)
        non_zero_count = np.count_nonzero(result)
        expected_count = np.pi * radius**2
        assert abs(non_zero_count - expected_count) < 100, "Radial cutoff area incorrect!"

    def test_no_radial_cutoff_none(self):
        """Test that None radial cutoff returns original mask"""
        test_mask = np.ones((100, 100))
        result = apply_radial_cutoff(test_mask, None)
        
        assert np.array_equal(result, test_mask), "Should return original!"

    def test_no_radial_cutoff_zero_radius(self):
        """Test that zero radius returns original mask"""
        test_mask = np.ones((100, 100))
        result = apply_radial_cutoff(test_mask, 0)
        
        assert np.array_equal(result, test_mask), "Zero radius should return original!"

    def test_radial_cutoff_negative_radius(self):
        """Test that negative radius returns original mask"""
        test_mask = np.ones((100, 100))
        result = apply_radial_cutoff(test_mask, -1)
        
        assert np.array_equal(result, test_mask), "Negative radius should return original!"

    def test_large_radius_cutoff(self):
        """Test large radius cutoff covering most of the image"""
        test_mask = np.ones((100, 100))
        result = apply_radial_cutoff(test_mask, 100)
        
        total_ones = np.count_nonzero(result)
        assert total_ones > 9000, "Large radius should cover most of image!"

    def test_radial_cutoff_circular_pattern(self):
        """Test that radial cutoff creates a circular pattern"""
        test_mask = np.ones((101, 101))  # Odd size for exact center
        radius = 25
        result = apply_radial_cutoff(test_mask, radius)
        
        center_y, center_x = 50, 50
        
        # Points within radius should be non-zero
        assert result[center_y, center_x] == 1.0
        assert result[center_y + 10, center_x] == 1.0
        assert result[center_y, center_x + 10] == 1.0
        
        # Points outside radius should be zero
        assert result[0, 0] == 0.0
        assert result[100, 100] == 0.0

    def test_radial_cutoff_with_nonbinary_mask(self):
        """Test radial cutoff on non-binary mask (should preserve multiplication)"""
        test_mask = np.full((100, 100), 5.0)
        radius = 30
        result = apply_radial_cutoff(test_mask, radius)
        
        center_y, center_x = 50, 50
        assert result[center_y, center_x] == 5.0, "Center should preserve value!"
        assert result[0, 0] == 0.0, "Corner should be zeroed!"


class TestIntegrationFOVAndRadialCutoff:
    """Integration tests for FOV and radial cutoff together"""

    def test_fov_crop_then_radial_cutoff(self):
        """Test applying FOV crop followed by radial cutoff"""
        # Create a test mask with varying values
        test_mask = np.arange(10000).reshape(100, 100).astype(float)
        
        # First apply FOV crop
        cropped_mask, fov = apply_fov_crop(test_mask, "25,25,75,75")
        assert cropped_mask.shape == (50, 50)
        
        # Then apply radial cutoff to the cropped region
        result = apply_radial_cutoff(cropped_mask, 20)
        assert result.shape == (50, 50)
        
        # Check that center is preserved and edges are zeroed
        center_y, center_x = 25, 25
        assert result[center_y, center_x] > 0, "Center should be non-zero!"
        assert result[0, 0] == 0.0, "Corner should be zero!"

    def test_radial_cutoff_then_optional_fov(self):
        """Test that operations can be done in different orders"""
        test_mask = np.ones((100, 100))
        
        # Apply radial cutoff first
        result1 = apply_radial_cutoff(test_mask, 30)
        
        # Then crop
        result2, _ = apply_fov_crop(result1, "20,20,80,80")
        
        # Should have circular region with square crop applied
        assert result2.shape == (60, 60)
        center_y, center_x = 30, 30
        # Center of cropped area should still be non-zero
        assert result2[center_y, center_x] == 1.0
