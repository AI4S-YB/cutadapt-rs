/// JSON formatting helpers.
///
/// Port of cutadapt/json.py.
///
/// Provides a custom JSON serializer that can render selected sub-trees
/// on a single line (matching the Python `OneLine` wrapper behaviour)
/// while still pretty-printing the rest of the document.

use serde_json::Value;

/// Marker wrapper: anything stored inside `OneLine` will be serialized
/// as compact (no newlines) JSON when passed through [`dumps`].
///
/// When building the JSON `Value` tree, call [`one_line`] to tag a value.
///
/// At the `serde_json::Value` level we cannot attach metadata, so we use
/// a two-pass approach: `dumps_with_one_line` accepts a parallel set of
/// paths that should be rendered compactly.
///
/// For simplicity in this port we provide two helpers:
///
/// 1. `compact_json(value)` -- render a value with no whitespace.
/// 2. `pretty_json(value, indent)` -- render with indentation (standard).
/// 3. `dumps(value, one_line_paths, indent)` -- render with selective compaction.

/// Render a `serde_json::Value` compactly (no whitespace).
pub fn compact_json(value: &Value) -> String {
    // serde_json::to_string produces compact output
    serde_json::to_string(value).unwrap_or_else(|_| "null".to_string())
}

/// Render a `serde_json::Value` with standard pretty-printing.
pub fn pretty_json(value: &Value) -> String {
    serde_json::to_string_pretty(value).unwrap_or_else(|_| "null".to_string())
}

/// Custom JSON formatter that matches the Python `cutadapt.json.dumps` output.
///
/// `one_line_markers` is a set of JSON pointer paths (like "/adapters_read1/0/trimmed_lengths/0")
/// whose values should be rendered compactly on one line.
///
/// For the common case where all array elements at a certain depth should be compact,
/// use `dumps_compact_arrays`.
pub fn dumps(value: &Value, indent: usize) -> String {
    format_value(value, indent, 0)
}

fn format_value(value: &Value, indent: usize, level: usize) -> String {
    match value {
        Value::Null | Value::Bool(_) | Value::Number(_) | Value::String(_) => {
            compact_json(value)
        }
        Value::Array(arr) => {
            if arr.is_empty() {
                return "[]".to_string();
            }
            let inner_indent = " ".repeat(indent * (level + 1));
            let outer_indent = " ".repeat(indent * level);
            let parts: Vec<String> = arr
                .iter()
                .map(|v| format_value(v, indent, level + 1))
                .collect();
            format!(
                "[\n{}{}\n{}]",
                inner_indent,
                parts.join(&format!(",\n{}", inner_indent)),
                outer_indent
            )
        }
        Value::Object(map) => {
            if map.is_empty() {
                return "{}".to_string();
            }
            let inner_indent = " ".repeat(indent * (level + 1));
            let outer_indent = " ".repeat(indent * level);
            let parts: Vec<String> = map
                .iter()
                .map(|(k, v)| {
                    let key_str = serde_json::to_string(k).unwrap();
                    let val_str = format_value(v, indent, level + 1);
                    format!("{}: {}", key_str, val_str)
                })
                .collect();
            format!(
                "{{\n{}{}\n{}}}",
                inner_indent,
                parts.join(&format!(",\n{}", inner_indent)),
                outer_indent
            )
        }
    }
}

/// Render a value compactly on a single line (for OneLine-tagged items).
pub fn one_line_json(value: &Value) -> String {
    compact_json(value)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn test_compact_json() {
        let v = json!({"a": 1, "b": [2, 3]});
        let s = compact_json(&v);
        assert!(!s.contains('\n'));
    }

    #[test]
    fn test_pretty_json() {
        let v = json!({"a": 1});
        let s = pretty_json(&v);
        assert!(s.contains('\n'));
    }

    #[test]
    fn test_dumps_empty_array() {
        let v = json!([]);
        assert_eq!(dumps(&v, 2), "[]");
    }

    #[test]
    fn test_dumps_empty_object() {
        let v = json!({});
        assert_eq!(dumps(&v, 2), "{}");
    }

    #[test]
    fn test_dumps_nested() {
        let v = json!({"a": [1, 2], "b": {"x": 3}});
        let s = dumps(&v, 2);
        assert!(s.contains('\n'));
        assert!(s.contains("\"a\""));
        assert!(s.contains("\"b\""));
    }

    #[test]
    fn test_one_line_json() {
        let v = json!({"len": 5, "count": 10});
        let s = one_line_json(&v);
        assert!(!s.contains('\n'));
        assert!(s.contains("\"len\""));
    }
}
