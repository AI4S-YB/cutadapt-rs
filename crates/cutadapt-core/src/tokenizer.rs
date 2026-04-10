/// Template string tokenization.
///
/// Port of cutadapt/tokenizer.py.

use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Token {
    String(String),
    Brace(String),
}

#[derive(Debug, Error)]
#[error("{0}")]
pub struct TokenizeError(String);

/// Tokenize a string by splitting on brace-delimited tokens.
///
/// # Examples
///
/// ```
/// use cutadapt_core::tokenizer::{tokenize_braces, Token};
///
/// let tokens: Vec<Token> = tokenize_braces("before {braced} after", '{', '}').unwrap();
/// assert_eq!(tokens, vec![
///     Token::String("before ".to_string()),
///     Token::Brace("braced".to_string()),
///     Token::String(" after".to_string()),
/// ]);
/// ```
pub fn tokenize_braces(s: &str, left: char, right: char) -> Result<Vec<Token>, TokenizeError> {
    if left == right {
        return Err(TokenizeError(
            "left and right must be unequal one-character strings".to_string(),
        ));
    }

    let mut tokens = Vec::new();
    let pattern = format!(
        "({}[^{}]*{})",
        regex::escape(&left.to_string()),
        regex::escape(&right.to_string()),
        regex::escape(&right.to_string())
    );
    let re = regex::Regex::new(&pattern).unwrap();

    let mut last_end = 0;
    for m in re.find_iter(s) {
        // Add any text before this match as a StringToken
        if m.start() > last_end {
            let before = &s[last_end..m.start()];
            if !before.is_empty() {
                if before.contains(left) {
                    return Err(TokenizeError(format!("Unexpected '{}' encountered", left)));
                }
                if before.contains(right) {
                    return Err(TokenizeError(format!(
                        "Unexpected '{}' encountered",
                        right
                    )));
                }
                tokens.push(Token::String(before.to_string()));
            }
        }

        // Extract the brace content (strip left and right delimiters)
        let matched = m.as_str();
        let inner = &matched[left.len_utf8()..matched.len() - right.len_utf8()];
        if inner.contains(left) {
            return Err(TokenizeError(format!("Unexpected '{}' encountered", left)));
        }
        if inner.contains(right) {
            return Err(TokenizeError(format!(
                "Unexpected '{}' encountered",
                right
            )));
        }
        tokens.push(Token::Brace(inner.to_string()));

        last_end = m.end();
    }

    // Add any remaining text after the last match
    if last_end < s.len() {
        let remaining = &s[last_end..];
        if !remaining.is_empty() {
            if remaining.contains(left) {
                return Err(TokenizeError(format!("Unexpected '{}' encountered", left)));
            }
            if remaining.contains(right) {
                return Err(TokenizeError(format!(
                    "Unexpected '{}' encountered",
                    right
                )));
            }
            tokens.push(Token::String(remaining.to_string()));
        }
    }

    Ok(tokens)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty() {
        let tokens = tokenize_braces("", '{', '}').unwrap();
        assert!(tokens.is_empty());
    }

    #[test]
    fn test_no_braces() {
        let tokens = tokenize_braces("hello world", '{', '}').unwrap();
        assert_eq!(tokens, vec![Token::String("hello world".to_string())]);
    }

    #[test]
    fn test_basic() {
        let tokens = tokenize_braces("before {braced} after", '{', '}').unwrap();
        assert_eq!(
            tokens,
            vec![
                Token::String("before ".to_string()),
                Token::Brace("braced".to_string()),
                Token::String(" after".to_string()),
            ]
        );
    }

    #[test]
    fn test_adjacent_braces() {
        let tokens = tokenize_braces("ab{cd}{ef}", '{', '}').unwrap();
        assert_eq!(
            tokens,
            vec![
                Token::String("ab".to_string()),
                Token::Brace("cd".to_string()),
                Token::Brace("ef".to_string()),
            ]
        );
    }

    #[test]
    fn test_custom_delimiters() {
        let tokens = tokenize_braces("ab(cd)ef", '(', ')').unwrap();
        assert_eq!(
            tokens,
            vec![
                Token::String("ab".to_string()),
                Token::Brace("cd".to_string()),
                Token::String("ef".to_string()),
            ]
        );
    }

    #[test]
    fn test_equal_delimiters() {
        let result = tokenize_braces("abc", '{', '{');
        assert!(result.is_err());
    }

    #[test]
    fn test_unmatched_left_brace() {
        let result = tokenize_braces("ab{cd", '{', '}');
        assert!(result.is_err());
    }

    #[test]
    fn test_unmatched_right_brace() {
        let result = tokenize_braces("ab}cd", '{', '}');
        assert!(result.is_err());
    }
}
