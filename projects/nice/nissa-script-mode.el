;;; nissa-script-mode.el --- Edit NICE scripts -*- lexical-binding: t; -*-

;; Completion reads the current buffer without executing its instructions

(require 'nissa-script-completion)
(require 'nissa-script-check)

(declare-function company-mode "company" (&optional arg))

(defvar nissa-script-mode-syntax-table
  (let ((table (make-syntax-table)))
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
  (let ((map (make-sparse-keymap)))
    (define-key map (kbd "TAB") #'completion-at-point)
    (define-key map (kbd "M-TAB") #'completion-at-point)
    (define-key map (kbd "C-c C-d") #'nissa-script-describe-entity)
    (define-key map (kbd "C-c C-v") #'flymake-start)
    map)
  "Keys for editing NICE scripts; TAB requests completion.")

;;;###autoload
(define-derived-mode nissa-script-mode prog-mode "NICE"
  "Edit NICE scripts with static completion from the NICE executable.
Set `nissa-script-completion-program' to the built executable.
\\<nissa-script-mode-map>Use \\[completion-at-point] to request completion.
The current buffer need not be saved.  Completion never executes its script."
  :syntax-table nissa-script-mode-syntax-table
  (setq-local font-lock-defaults '(nissa-script-mode-font-lock-keywords))
  (setq-local comment-start "// ")
  (setq-local comment-end "")
  (setq-local comment-start-skip "//+\\s-*\\|/\\*+\\s-*")
  (nissa-script-completion-mode 1)
  (nissa-script-check-mode 1)
  (when (require 'company nil t)
    (setq-local company-backends '(company-capf))
    (setq-local company-minimum-prefix-length 1)
    (setq-local company-idle-delay 0.2)
    (company-mode 1)))

(provide 'nissa-script-mode)
;;; nissa-script-mode.el ends here
