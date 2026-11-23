#ifndef _NICE_SCRIPT_SOURCE_HPP
#define _NICE_SCRIPT_SOURCE_HPP

#include <algorithm>
#include <iterator>
#include <memory>
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
  
  /// Source line owned by the AST after the input buffer is released
  std::string lineText{};
  
  /// Token width for the runtime error underline
  size_t extent{1};
  
  /// Copy a parser location and its source excerpt into AST storage
  static ScriptLocation fromPosition(const pp::internal::SourceLocation& position,
                                     const size_t extent)
  {
    return {position.source.empty()?"<script>":std::string(position.source),
            position.cursor,std::string(position.lineText),extent};
  }
  
  /// Format a location suitable for editor navigation and diagnostics
  std::string describe() const
  {
    return source+":"+std::to_string(cursor.line)+":"+std::to_string(cursor.column);
  }
  
  /// Report an evaluation error at this script position
  template <typename...Args>
  [[noreturn]] void error(Args&&...args) const
  {
    pp::internal::SourceLocation{source,cursor,lineText,extent}.error(std::forward<Args>(args)...);
  }
};

/// Supply runtime errors in nested operators and host functions with their script position
struct ScriptDiagnosticContext
{
  /// Borrow the AST's source label for this evaluation scope
  pp::internal::SourceLocation position;
  
  /// Restore the enclosing expression's diagnostic position on scope exit
  pp::internal::ScopedErrorLocation scope;
  
  /// Nodes without a position inherit the enclosing expression's context
  explicit ScriptDiagnosticContext(const ScriptLocation* location) :
    position{location?location->source:std::string_view{},location?location->cursor:pp::internal::SourceCursor{},
             location?location->lineText:std::string_view{},location?location->extent:1},
    scope(location?&position:pp::internal::errorLocation)
  {
  }
};

/// Owned invocation location retained when a reference outlives its script call
struct ScriptCallSite
{
  /// Callable name as seen by its caller
  std::string name;
  
  /// Owned invocation position and source excerpt
  ScriptLocation position;
};

/// Linked script invocation frames, independent of the C++ backtrace
struct ScriptCallContext
{
  /// Innermost invocation active on the current interpreter thread
  inline static thread_local const ScriptCallContext* active=nullptr;
  
  /// Saved invocation chain selected by a deferred reference read or write
  inline static thread_local const std::vector<ScriptCallSite>* retained=nullptr;
  
  /// Enclosing script invocation
  const ScriptCallContext* previous;
  
  /// Restore an embedding application's diagnostic footer on return
  void (*previousFooter)(std::ostream&);
  
  /// Script-visible callable name, or an anonymous-call placeholder
  std::string_view name;
  
  /// Invocation position, borrowing the enclosing AST's retained strings
  pp::internal::SourceLocation position;
  
  /// Register this invocation before evaluating arguments and defaults
  explicit ScriptCallContext(const std::string_view name) :
    previous(active),
    previousFooter(pp::internal::errorFooter),
    name(name.empty()?"<callable>":name),
    position(pp::internal::errorLocation?*pp::internal::errorLocation:pp::internal::SourceLocation{})
  {
    active=this;
    pp::internal::errorFooter=write;
  }
  
  /// Restore the caller without leaving frames from completed invocations
  ~ScriptCallContext()
  {
    active=previous;
    pp::internal::errorFooter=previousFooter;
  }
  
  /// Each invocation owns exactly one stack entry
  ScriptCallContext(const ScriptCallContext&)=delete;
  
  /// An active entry cannot be reassigned
  ScriptCallContext& operator=(const ScriptCallContext&)=delete;
  
  /// Print innermost-to-outermost script calls only when an error is emitted
  static void write(std::ostream& output)
  {
    output<<"Script call stack (innermost first):\n";
    if(retained)
      for(const auto& frame : *retained)
        output<<"  "<<frame.name<<" at "<<frame.position.describe()<<'\n';
    else
      for(auto* frame=active;frame;frame=frame->previous)
        {
          output<<"  "<<frame->name<<" at ";
          if(frame->position.source.empty())
            output<<"<script>";
          else
            output<<frame->position.source;
          output<<':'<<frame->position.cursor.line<<':'<<frame->position.cursor.column<<'\n';
        }
  }
};

/// Restore a returned reference's invocation chain only while it is accessed
struct ScriptDeferredCallContext
{
  /// Restore any surrounding deferred access
  const std::vector<ScriptCallSite>* previous;
  
  /// Restore the enclosing diagnostic footer
  void (*previousFooter)(std::ostream&);
  
  /// Select a saved chain without changing normal calls when none was retained
  explicit ScriptDeferredCallContext(const std::vector<ScriptCallSite>* calls) :
    previous(ScriptCallContext::retained),
    previousFooter(pp::internal::errorFooter)
  {
    if(calls)
      {
        ScriptCallContext::retained=calls;
        pp::internal::errorFooter=ScriptCallContext::write;
      }
  }
  
  /// Leave no saved diagnostic chain after the access completes
  ~ScriptDeferredCallContext()
  {
    ScriptCallContext::retained=previous;
    pp::internal::errorFooter=previousFooter;
  }
  
  /// A diagnostic scope has one owner
  ScriptDeferredCallContext(const ScriptDeferredCallContext&)=delete;
  
  /// An active diagnostic scope cannot be reassigned
  ScriptDeferredCallContext& operator=(const ScriptDeferredCallContext&)=delete;
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
        
        const size_t lineBegin=piece.begin+cursor.offset-(cursor.column-1);
        size_t lineEnd=lineBegin;
        while(lineEnd<piece.begin+piece.size and text[lineEnd]!='\n' and text[lineEnd]!='\r')
          lineEnd++;
        result={piece.origin,cursor,std::string_view(text).substr(lineBegin,lineEnd-lineBegin)};
      }
    
    return result;
  }
};

#endif
