#ifndef _NICE_SCRIPT_SOURCE_HPP
#define _NICE_SCRIPT_SOURCE_HPP

#include <algorithm>
#include <iterator>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <parsePact.hpp>

/// Owned declaration position, independent of parser and source-text lifetimes
struct ScriptLocation
{
  /// File path, argument label, or <script> for an unnamed input
  std::string source;
  
  /// Zero-based byte offset and one-based line and byte column in that source
  pp::internal::SourceCursor cursor;
  
  /// Format a location suitable for editor navigation and diagnostics
  std::string describe() const
  {
    return source+":"+std::to_string(cursor.line)+":"+std::to_string(cursor.column);
  }
  
  /// Report an evaluation error at this script position
  template <typename...Args>
  [[noreturn]] void error(Args&&...args) const
  {
    pp::internal::errorEmitter(describe(),": ",std::forward<Args>(args)...);
  }
};

/// Joined program text with original-source positions preserved across pieces
struct ScriptSource
{
  /// One contiguous source piece and its precomputed line starts
  struct Piece
  {
    /// Beginning of this piece in the joined text
    size_t begin;
    
    /// Original byte count, excluding an inserted terminating newline
    size_t size;
    
    /// Owned original source label
    std::string origin;
    
    /// Local cursors just after newlines, including the beginning of the piece
    std::vector<pp::internal::SourceCursor> lines;
  };
  
  /// Joined input retained throughout parsing and parse-tree execution
  std::string text;
  
  /// Source pieces in increasing joined-byte order
  std::vector<Piece> pieces;
  
  /// Append a piece and terminate its trailing line comment when necessary
  void append(const std::string_view input,
              const std::string& origin)
  {
    Piece piece{.begin=text.size(),.size=input.size(),.origin=origin,.lines={{}}};
    text+=input;
    
    if(input.empty() or (input.back()!='\n' and input.back()!='\r'))
      text+='\n';
    
    pp::internal::SourceCursor cursor;
    for(size_t i=piece.begin;i<text.size();i++)
      {
        cursor.advance(text[i]);
        
        if(text[i]=='\r' or text[i]=='\n')
          {
            if(piece.lines.back().line==cursor.line)
              piece.lines.back()=cursor;
            else
              piece.lines.push_back(cursor);
          }
      }
    
    pieces.push_back(std::move(piece));
  }
  
  /// Resolve a joined cursor to a borrowed label and local cursor
  /// Returned labels remain valid until this source is modified or destroyed
  pp::internal::SourceLocation locate(const pp::internal::SourceCursor& joined) const
  {
    pp::internal::SourceLocation result{{},joined};
    
    if(not pieces.empty())
      {
        const auto next=std::upper_bound(pieces.begin(),pieces.end(),joined.offset,
          [](const size_t offset,
             const Piece& piece)
          {
            return offset<piece.begin;
          });
        const auto& piece=*(next==pieces.begin()?next:std::prev(next));
        const size_t offset=std::min(std::min(joined.offset,text.size())-piece.begin,piece.size);
        const auto line=std::upper_bound(piece.lines.begin(),piece.lines.end(),offset,
          [](const size_t position,
             const pp::internal::SourceCursor& cursor)
          {
            return position<cursor.offset;
          });
        auto cursor=*std::prev(line);
        while(cursor.offset<offset)
          cursor.advance(text[piece.begin+cursor.offset]);
        
        result={piece.origin,cursor};
      }
    
    return result;
  }
};

#endif
