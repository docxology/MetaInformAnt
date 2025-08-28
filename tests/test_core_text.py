from __future__ import annotations

import pytest
from metainformant.core import text as core_text


class TestSlugify:
    """Comprehensive tests for text slugification."""
    
    def test_slugify_basic(self) -> None:
        """Test basic slugification functionality."""
        assert core_text.slugify("Hello, World!") == "hello-world"
        assert core_text.slugify("  A  B  C  ") == "a-b-c"

    def test_slugify_edge_cases(self) -> None:
        """Test edge cases for slugification."""
        # Empty string
        assert core_text.slugify("") == ""
        
        # Only whitespace
        assert core_text.slugify("   ") == ""
        
        # Only punctuation
        assert core_text.slugify("!@#$%^&*()") == ""
        
        # Mixed case with numbers
        assert core_text.slugify("Test123 File") == "test123-file"
        
        # Unicode characters (accents are removed by the regex)
        result = core_text.slugify("Café résumé naïve")
        assert result == "caf-rsum-nave"  # Accented chars are removed
        
        # Multiple consecutive spaces/punctuation
        assert core_text.slugify("a...b___c---d") == "abc-d"  # Keeps last dash
        
        # Leading/trailing separators
        assert core_text.slugify("---test---") == "test"

    def test_slugify_very_long_input(self) -> None:
        """Test slugification with very long input."""
        long_string = "word " * 1000
        result = core_text.slugify(long_string)
        # Should handle long strings without issue
        assert len(result) > 0
        assert result.startswith("word")
        assert not result.endswith("-")

    def test_slugify_special_characters(self) -> None:
        """Test handling of various special characters."""
        # The actual slugify implementation removes all non-alphanumeric chars except dashes
        test_cases = [
            ("file (1).txt", "file-1txt"),  # Periods are removed
            ("test@domain.com", "testdomaincom"),  # @ and . are removed
            ("100% complete", "100-complete"),
            ("C++ programming", "c-programming"),  # + signs are removed
            ("file [draft].doc", "file-draftdoc")  # Brackets and period removed
        ]
        
        for input_str, expected in test_cases:
            assert core_text.slugify(input_str) == expected

    def test_slugify_non_string_input(self) -> None:
        """Test error handling for non-string input."""
        with pytest.raises((TypeError, AttributeError)):
            core_text.slugify(123)
            
        with pytest.raises((TypeError, AttributeError)):
            core_text.slugify(None)


class TestNormalizeWhitespace:
    """Comprehensive tests for whitespace normalization."""
    
    def test_normalize_whitespace_basic(self) -> None:
        """Test basic whitespace normalization."""
        s = "a\t\t b\n\n c"
        assert core_text.normalize_whitespace(s) == "a b c"

    def test_normalize_whitespace_edge_cases(self) -> None:
        """Test edge cases for whitespace normalization."""
        # Empty string
        assert core_text.normalize_whitespace("") == ""
        
        # Only whitespace
        assert core_text.normalize_whitespace("   \t\n\r  ") == ""
        
        # Single word
        assert core_text.normalize_whitespace("hello") == "hello"
        
        # Mixed whitespace types
        assert core_text.normalize_whitespace("a\n\r\t b\v\f c") == "a b c"
        
        # Leading/trailing whitespace
        assert core_text.normalize_whitespace("  hello world  ") == "hello world"

    def test_normalize_whitespace_unicode(self) -> None:
        """Test whitespace normalization with unicode characters."""
        # Unicode spaces (non-breaking space, em space, etc.)
        unicode_text = "hello\u00A0world\u2003test"
        result = core_text.normalize_whitespace(unicode_text)
        # Should normalize various unicode whitespace
        assert "hello" in result and "world" in result and "test" in result

    def test_normalize_whitespace_very_long(self) -> None:
        """Test normalization with very long strings."""
        long_string = "a" + " " * 1000 + "b" + "\t" * 500 + "c"
        result = core_text.normalize_whitespace(long_string)
        assert result == "a b c"

    def test_normalize_whitespace_preserve_structure(self) -> None:
        """Test that normalization preserves word order."""
        text = "  first   second\t\tthird\n\nfourth  "
        result = core_text.normalize_whitespace(text)
        assert result == "first second third fourth"


class TestSafeFilename:
    """Comprehensive tests for safe filename generation."""
    
    def test_safe_filename_basic(self) -> None:
        """Test basic safe filename functionality."""
        name = "My*Weird:File?.txt"
        safe = core_text.safe_filename(name)
        assert all(ch not in safe for ch in ['*', ':', '?'])
        assert safe.endswith(".txt")

    def test_safe_filename_dangerous_characters(self) -> None:
        """Test removal of dangerous filesystem characters."""
        dangerous = r'file<>:"/\|?*name.txt'
        safe = core_text.safe_filename(dangerous)
        
        # Should remove all dangerous characters
        dangerous_chars = r'<>:"/\|?*'
        for char in dangerous_chars:
            assert char not in safe
            
        # The result should be "name.txt" since dangerous chars are removed by slugify
        assert safe == "name.txt"
        assert "name" in safe
        assert ".txt" in safe

    def test_safe_filename_edge_cases(self) -> None:
        """Test edge cases for safe filename generation."""
        # Empty string
        assert core_text.safe_filename("") == ""
        
        # Only dangerous characters
        dangerous_only = r'<>:"/\|?*'
        result = core_text.safe_filename(dangerous_only)
        # Should return empty string since all chars are stripped
        assert result == ""
        
        # Very long filename - the implementation doesn't truncate
        long_name = "very_long_filename_" * 20 + ".txt"
        safe_long = core_text.safe_filename(long_name)
        # Document that it doesn't truncate (implementation doesn't enforce limits)
        assert safe_long.endswith(".txt")
        # Underscores are removed, not converted to dashes
        assert "verylongfilename" in safe_long

    def test_safe_filename_preserve_extension(self) -> None:
        """Test that file extensions are preserved correctly."""
        test_cases = [
            ("document.pdf", ".pdf"),
            ("image.jpg", ".jpg"),
            ("archive.tar.gz", ".gz"),  # Should preserve last extension
            ("no_extension", ""),
            ("multiple.dots.in.name.txt", ".txt")
        ]
        
        for filename, expected_ext in test_cases:
            safe = core_text.safe_filename(filename)
            if expected_ext:
                assert safe.endswith(expected_ext)

    def test_safe_filename_reserved_names(self) -> None:
        """Test handling of reserved filename patterns."""
        # Windows reserved names - the current implementation doesn't handle these specially
        # It just applies slugify which lowercases them
        reserved_names = ["CON", "PRN", "AUX", "NUL", "COM1", "LPT1"]
        
        for name in reserved_names:
            safe = core_text.safe_filename(name + ".txt")
            # The implementation currently just lowercases, doesn't handle reserved names
            assert safe == (name + ".txt").lower()
            # This documents current behavior - could be enhanced in future

    def test_safe_filename_unicode(self) -> None:
        """Test handling of unicode characters in filenames."""
        unicode_name = "café_résumé_naïve.txt"
        safe = core_text.safe_filename(unicode_name)
        
        # Should handle unicode gracefully (behavior may vary by implementation)
        assert len(safe) > 0
        assert ".txt" in safe

    def test_safe_filename_whitespace(self) -> None:
        """Test handling of whitespace in filenames."""
        name_with_spaces = "  file with spaces  .txt"
        safe = core_text.safe_filename(name_with_spaces)
        
        # Should handle leading/trailing spaces
        assert not safe.startswith(" ")
        assert not safe.endswith(" .txt")
        assert "file" in safe and "with" in safe and "spaces" in safe


class TestTextUtilsIntegration:
    """Integration tests for text utility functions."""
    
    def test_text_processing_pipeline(self) -> None:
        """Test combining multiple text processing functions."""
        original = "  My*Weird<File>Name?.doc  "
        
        # Step 1: Normalize whitespace
        normalized = core_text.normalize_whitespace(original)
        assert normalized == "My*Weird<File>Name?.doc"
        
        # Step 2: Make safe filename
        safe = core_text.safe_filename(normalized)
        dangerous_chars = r'*<>?'
        for char in dangerous_chars:
            assert char not in safe
        assert ".doc" in safe
        
        # Step 3: Create slug
        slug = core_text.slugify(normalized.replace(".doc", ""))
        assert slug == "myweirdfilename"  # All special chars are removed

    def test_empty_input_consistency(self) -> None:
        """Test that all functions handle empty input consistently."""
        empty_inputs = ["", "   ", "\t\n\r"]
        
        for empty in empty_inputs:
            slug = core_text.slugify(empty)
            normalized = core_text.normalize_whitespace(empty)
            safe = core_text.safe_filename(empty)
            
            # All should return empty string for empty input
            assert slug == ""
            assert normalized == ""
            assert safe == ""

    def test_round_trip_compatibility(self) -> None:
        """Test that functions work well together in various orders."""
        test_string = "Test File (Version 1.2).txt"
        
        # Different processing orders should be stable
        order1 = core_text.safe_filename(core_text.normalize_whitespace(test_string))
        order2 = core_text.normalize_whitespace(core_text.safe_filename(test_string))
        
        # Both should be safe filenames
        assert ".txt" in order1
        assert ".txt" in order2
        # Both should preserve the essential content (lowercased by slugify)
        assert "test" in order1 and "file" in order1  # Lowercase due to slugify
        assert "test" in order2 and "file" in order2


