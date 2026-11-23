;;; nissa-script-mode.el --- Edit NICE scripts -*- lexical-binding: t; -*-

;; Completion reads the current buffer without executing its instructions

(require 'nissa-script-completion)
(require 'nissa-script-check)
(require 'cc-mode)

(declare-function company-mode "company" (&optional arg))

(defvar nissa-script-mode-syntax-table
  (let ((table (copy-syntax-table c-mode-syntax-table)))
    (modify-syntax-entry ?_ "w" table)
    (modify-syntax-entry ?/ ". 124b" table)
    (modify-syntax-entry ?* ". 23" table)
    (modify-syntax-entry ?\n "> b" table)
    (modify-syntax-entry ?\r "> b" table)
    (modify-syntax-entry ?\" "\"" table)
    (modify-syntax-entry ?\\ "\\" table)
    (modify-syntax-entry ?\( "()" table)
    (modify-syntax-entry ?\) ")(" table)
    (modify-syntax-entry ?\[ "(]" table)
    (modify-syntax-entry ?\] ")[" table)
    (modify-syntax-entry ?{ "(}" table)
    (modify-syntax-entry ?} "){" table)
    table)
  "Syntax of NICE identifiers, strings, comments, and delimiters.")

(defvar nissa-script-mode-font-lock-keywords
  `((,(regexp-opt '("fun" "struct" "lambda" "if" "else" "for" "return") 'symbols)
     . font-lock-keyword-face)
    ("\\_<\\(?:fun\\|struct\\)\\_>\\s-+\\([[:word:]_]+\\)"
     1 font-lock-function-name-face))
  "Basic highlighting for NICE declarations and control flow.")

(defvar nissa-script-mode-map
  (make-sparse-keymap)
  "Keys for editing NICE scripts; TAB indents and then offers completion.")

;; Refresh bindings when this file is reloaded into an existing Emacs session
(define-key nissa-script-mode-map (kbd "TAB") #'indent-for-tab-command)
(define-key nissa-script-mode-map (kbd "M-TAB") #'completion-at-point)
(define-key nissa-script-mode-map (kbd "C-c C-d") #'nissa-script-describe-entity)
(define-key nissa-script-mode-map (kbd "C-c C-v") #'flymake-start)

;;;###autoload
(define-derived-mode nissa-script-mode c-mode "NICE"
  "Edit NICE scripts with C indentation and static completion.
Indentation follows the C style selected by `c-default-style'; customize
`c-basic-offset' or `c-offsets-alist' in `nissa-script-mode-hook' as needed.
Set `nissa-script-completion-program' to the built executable.
TAB indents first and offers completion when indentation is already correct.
\\<nissa-script-mode-map>Use \\[completion-at-point] for completion directly.
The current buffer need not be saved.  Completion never executes its script."
  :syntax-table nissa-script-mode-syntax-table
  ;; Parent hooks run after this body.  Keep their editing customizations, but
  ;; do not start C language servers for NICE's separate script language
  (setq-local c-mode-hook (copy-sequence c-mode-hook))
  (remove-hook 'c-mode-hook #'lsp t)
  (remove-hook 'c-mode-hook #'lsp-deferred t)
  ;; NICE structs end at their closing brace, with no C declarator or semicolon
  ;; afterward; otherwise CC Mode treats the following definition as a continuation
  (setq-local c-opt-block-decls-with-vars-key nil)
  (setq-local tab-always-indent 'complete)
  (setq-local font-lock-defaults '(nissa-script-mode-font-lock-keywords))
  (setq-local comment-start "// ")
  (setq-local comment-end "")
  (setq-local comment-start-skip "//+\\s-*\\|/\\*+\\s-*")
  (nissa-script-completion-mode 1)
  ;; C indentation does not make a NICE buffer a C compilation unit
  (remove-hook 'flymake-diagnostic-functions #'flymake-cc t)
  (nissa-script-check-mode 1)
  (when (require 'company nil t)
    (setq-local company-backends '(company-capf))
    (setq-local company-minimum-prefix-length 1)
    (setq-local company-idle-delay 0.2)
    (company-mode 1)))

(provide 'nissa-script-mode)
;;; nissa-script-mode.el ends here
