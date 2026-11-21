R"NASTY_HELP(nASTy script language
=====================

Programs and values
  A program is a sequence of statements; outer braces are unnecessary.
  End expressions and assignments with a semicolon. Braces create an inner scope.
  Extra semicolons are empty statements and have no effect.
  Assignments create variables: count=3; step=0.1; label="example";
  Values include integers, floating-point numbers, strings, lists, functions,
  script structs, and application-defined host objects.
  Use // or /* ... */ for comments. Immediately preceding /// or //! lines
  document a function, struct, or assigned variable.

Expressions
  Arithmetic: +, -, *, /, % and unary + or -; parentheses group expressions.
  Comparison: <, <=, >, >=, ==, !=. Logical: !, && (and), || (or).
  Logical and/or short-circuit: the right operand is evaluated only if needed.
  Multiplication/division/remainder bind before addition/subtraction; comparison
  binds before logical and, which binds before logical or. Assignment binds last.
  Arithmetic operators at the same precedence associate to the left.
  Assignments: =, +=, -=, *=, /=. Postfix updates: count++; count--;
  Index a list with values[index] (zero-based); access a member with object.member.
  Host objects can supply application-specific operators and member access.

Calls and lists
  Call functions with positional or named arguments: f(3,.scale=2);
  Named arguments use .name=value; omitted optional parameters use their defaults.
  makeVec(1,2,3) creates a list. seq(max), seq(min,max), and seq(min,step,max)
  create numeric sequences. Use help("seq") for its registered description.
  print("value=",count,"\n"); prints without an automatic final newline.

Control flow
  if(condition) { statements } else { statements }
  for(i=0;i<3;i++) { statements }
  return expression; returns from a function.

Functions and structs
  fun scale(value,factor=2) { return value*factor; }
  twice=lambda(value) { return 2*value; };
  Prefix a parameter with & to accept a reference: fun update(&value) { value++; }
  A final ... parameter collects additional arguments in __VA__ARGS__; &... collects
  references instead. Use &value when explicitly passing a reference expression.
  struct Pair { first; second=0; };
  pair=Pair(1,.second=2); pair.first=3;

Discovering the available interface
  help();                  Show this language overview
  help("name");            Show the nearest visible entity's signature and description
  listFunctions();          List visible functions and constructors
  listVariables();          List visible variables and host objects
  listEntities();           List all visible entities
  Listings default to .verbose=1; .verbose=0 prints names only.
  In signatures, [...] marks optional arguments, & marks references, and ... marks
  variadic arguments. Defaults appear inline when available; <expression> means
  a default is not evaluated for documentation.
  describe("name","text"); attaches documentation to an existing visible binding.
  Inner bindings hide outer bindings. Application-specific functions, operators,
  and constructors are documented through their registered descriptions.
)NASTY_HELP"
