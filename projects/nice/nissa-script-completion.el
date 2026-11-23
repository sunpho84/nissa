;;; nissa-script-completion.el --- CAPF adapter for nASTy queries -*- lexical-binding: t; -*-

;; Backend-independent adapter for formatScriptCompletionForEmacs
;; Querying must not execute the input script

(require 'subr-x)
(require 'cl-lib)
(require 'eldoc)

(declare-function yas-expand-snippet "yasnippet" (snippet &optional start end expand-env))
(declare-function yas-minor-mode "yasnippet" (&optional arg))

(defgroup nissa-script nil
  "Editing NICE scripts."
  :group 'languages)

(defcustom nissa-script-completion-program "NICE"
  "NICE executable used for completion queries."
  :type 'string
  :group 'nissa-script)

(require 'nissa-script-session)

(defcustom nissa-script-completion-context nil
  "Ordered program pieces preceding the current buffer.
Each entry is (file . PATH) or (script . TEXT).  Paths are relative to the
script buffer's `default-directory'.  Files are read for static completion
only; no script is executed.  Set this buffer-locally for a particular run."
  :type '(repeat (cons (choice (const file) (const script)) string))
  :group 'nissa-script)

(defun nissa-script-completion-context-arguments ()
  "Return command-line arguments for the configured preceding program pieces."
  (cl-mapcan
                  (lambda (piece)
                    (list (pcase (car piece)
                            ('file "--context-file")
                            ('script "--context-script")
                            (_ (error "Unknown NICE completion context kind: %S" (car piece))))
                          (cdr piece)))
                  nissa-script-completion-context))

(defun nissa-script-completion-query (text byte-offset)
  "Complete TEXT at BYTE-OFFSET, reusing this buffer's worker by default."
  (condition-case error-data
      (if nissa-script-use-persistent-session
          (nissa-script-session-query text byte-offset)
        (nissa-script-session-stop)
        (nissa-script-completion-query-once text byte-offset))
    (error (message "NICE completion failed: %s" (error-message-string error-data)) nil)))

(defun nissa-script-completion-query-once (text byte-offset)
  "Send TEXT and BYTE-OFFSET to NICE without evaluating the script."
  (let ((program nissa-script-completion-program)
        (context (nissa-script-completion-context-arguments))
        (origin (or buffer-file-name "<script>"))
        (directory default-directory)
        (errors (make-temp-file "nissa-completion-"))
        (coding-system-for-read 'utf-8-unix)
        (coding-system-for-write 'utf-8-unix))
    (unwind-protect
        (with-temp-buffer
          (insert text)
          (let* ((default-directory directory)
                 (status (apply #'call-process-region
                                (point-min) (point-max) program
                                t (list t errors) nil
                                "--complete-at" (number-to-string byte-offset)
                                "--completion-format" "emacs"
                                "--completion-source" origin context)))
            (if (equal status 0)
                (progn (goto-char (point-min)) (read (current-buffer)))
              (message "NICE completion failed: %s"
                       (with-temp-buffer
                         (insert-file-contents errors)
                         (string-trim (buffer-string))))
              nil)))
      (delete-file errors))))

(defvar-local nissa-script-completion-query-function #'nissa-script-completion-query
  "Function accepting source TEXT and a zero-based UTF-8 BYTE-OFFSET.
Return nil, or (BEGIN END CANDIDATES), where BEGIN and END are zero-based
UTF-8 byte offsets and each candidate is
(NAME KIND DOCUMENTATION SOURCE LINE COLUMN OFFSET).
The default provider shares a persistent worker with Flymake, before normal
NISSA/MPI initialization, and never executes TEXT.
If using a textual response, read its Lisp data with `read', never `eval'.")

(defun nissa-script-completion-signature (entry)
  "Return ENTRY's callable signature, or nil when none is available."
  (let* ((kind (nth 1 entry))
         (first-line (string-trim
                      (or (car (split-string (or (nth 2 entry) "") "\n" t)) ""))))
    (when (and (string-match-p "\\(function\\|constructor\\|lambda\\)" kind)
               (string-suffix-p ")" first-line)
               (string-match "\\`[[:word:]_]+[[:space:]]*(" first-line))
      (concat (car entry) (substring first-line (1- (match-end 0)))))))

(defun nissa-script-completion-doc-buffer (entry)
  "Return a documentation buffer for ENTRY without displaying it."
  (when-let* ((documentation (nth 2 entry))
              ((not (string-empty-p documentation))))
    (with-current-buffer (get-buffer-create "*NICE completion documentation*")
      (let ((inhibit-read-only t))
        (erase-buffer)
        (insert documentation)
        (goto-char (point-min))
        (special-mode))
      (current-buffer))))

(defun nissa-script-completion-summary (entry)
  "Return ENTRY's signature and a short description for the echo area."
  (let* ((documentation (or (nth 2 entry) ""))
         (callable (nissa-script-completion-signature entry))
         (signature (or callable (car entry)))
         (first-line (car (split-string documentation "\n" t)))
         (description
          (if (string-match "Description:\n" documentation)
              (substring documentation (match-end 0))
            (mapconcat #'identity
                       (cl-remove-if
                        (lambda (line)
                          (or (and callable (equal line first-line))
                              (string-match-p "^[[:space:]]*\\(?:Kind\\|Location\\):" line)))
                        (split-string documentation "\n" t))
                       " "))))
    (concat signature
            (if (string-empty-p (string-trim description)) ""
              (concat " — " (string-join (split-string description) " "))))))

(defun nissa-script-completion-required-arguments (entry)
  "Return required parameter names from ENTRY's signature.
Optional and variadic parameters are omitted from inserted calls."
  (when-let* ((signature (nissa-script-completion-signature entry))
              (open (string-match "(" signature)))
    (let ((begin (1+ open)) (depth 0) (quoted nil) (escaped nil) arguments)
      (cl-loop for i from begin below (length signature)
               for c = (aref signature i)
               do (cond
                   (quoted
                    (when (and (= c ?\") (not escaped)) (setq quoted nil))
                    (setq escaped (and (= c ?\\) (not escaped))))
                   ((= c ?\") (setq quoted t))
                   ((memq c '(?\[ ?\()) (setq depth (1+ depth)))
                   ((and (memq c '(?\] ?\))) (> depth 0)) (setq depth (1- depth)))
                   ((and (= depth 0) (memq c '(?, ?\))))
                    (let ((argument (string-trim (substring signature begin i))))
                      (when (string-match "\\`&?\\([[:alpha:]_][[:alnum:]_]*\\)\\'" argument)
                        (push (match-string 1 argument) arguments)))
                    (setq begin (1+ i)))))
      (nreverse arguments))))

(defun nissa-script-completion-insert-call (entry)
  "Insert ENTRY's call template unless an argument list already follows point."
  (when (and (nissa-script-completion-signature entry)
             (not (save-excursion
                    (skip-chars-forward " \t\n")
                    (looking-at-p "("))))
    (let ((arguments (nissa-script-completion-required-arguments entry)))
      (if (and arguments (require 'yasnippet nil t))
          (let ((index 0))
            (yas-minor-mode 1)
            (yas-expand-snippet
             (concat "("
                     (mapconcat (lambda (name)
                                  (format "${%d:%s}" (cl-incf index) name))
                                arguments ", ")
                     ")$0")))
        (insert "(" (string-join arguments ", ") ")")
        (if arguments
            (progn
              (backward-char (1+ (length (string-join arguments ", "))))
              (push-mark (+ (point) (length (car arguments))) t t))
          (backward-char 1))))))

(defvar-local nissa-script-completion-entity-cache nil
  "Last entity documentation query and its result.")

(defun nissa-script-completion-entity (bounds)
  "Query the entity occupying BOUNDS without changing point."
  (when (and bounds nissa-script-completion-query-function)
    (save-restriction
      (widen)
      (let* ((name (buffer-substring-no-properties (car bounds) (cdr bounds)))
             (key (list (buffer-chars-modified-tick) (cdr bounds) buffer-file-name
                        nissa-script-completion-context nissa-script-completion-program
                        nissa-script-completion-query-function default-directory)))
        (unless (and (not (and nissa-script-use-persistent-session
                              (eq nissa-script-completion-query-function #'nissa-script-completion-query)))
                     (equal key (car nissa-script-completion-entity-cache)))
          (let* ((text (buffer-substring-no-properties (point-min) (point-max)))
                 (offset (length (encode-coding-string
                                  (buffer-substring-no-properties (point-min) (cdr bounds))
                                  'utf-8 t)))
                 (reply (funcall nissa-script-completion-query-function text offset)))
            (setq nissa-script-completion-entity-cache
                  (list key (assoc-string name (nth 2 reply))))))
        (cadr nissa-script-completion-entity-cache)))))

(defun nissa-script-completion-entity-at-point ()
  "Find metadata for the entity at point, or for the enclosing call."
  (save-excursion
    (save-match-data
      (unless (nth 8 (syntax-ppss))
        (or (nissa-script-completion-entity (bounds-of-thing-at-point 'symbol))
            (when-let* ((open (nth 1 (syntax-ppss))))
              (goto-char open)
              (when (looking-at-p "(")
                (skip-chars-backward " \t\n")
                (nissa-script-completion-entity (bounds-of-thing-at-point 'symbol)))))))))

(defun nissa-script-completion-eldoc (callback &rest _ignored)
  "Provide CALLBACK with documentation for the current NICE entity."
  (when-let* ((entry (nissa-script-completion-entity-at-point)))
    (funcall callback (or (nth 2 entry) (car entry))
             :thing (car entry) :echo (nissa-script-completion-summary entry))))

(defun nissa-script-describe-entity ()
  "Display full documentation for the NICE entity at point."
  (interactive)
  (if-let* ((entry (nissa-script-completion-entity-at-point))
             (buffer (nissa-script-completion-doc-buffer entry)))
      (display-buffer buffer)
    (message "No NICE entity documentation at point")))

(defun nissa-script-completion-at-point ()
  "Complete nASTy identifiers using the configured non-executing query.
Install buffer-locally in `completion-at-point-functions'."
  (when nissa-script-completion-query-function
    (save-restriction
      (widen)
      (let* ((text (buffer-substring-no-properties (point-min) (point-max)))
             (offset (length (encode-coding-string
                              (buffer-substring-no-properties (point-min) (point))
                              'utf-8 t)))
             (reply (funcall nissa-script-completion-query-function text offset)))
        (when reply
          (pcase-let* ((`(,begin ,end ,entries) reply)
                       (size (length (encode-coding-string text 'utf-8 t))))
            (when (and (integerp begin) (integerp end)
                       (<= 0 begin offset end size))
              ;; The external protocol always uses UTF-8, independently of
              ;; Emacs's internal byte representation
              (let* ((bytes (encode-coding-string text 'utf-8 t))
                     (beg-pos (1+ (length (decode-coding-string
                                          (substring bytes 0 begin) 'utf-8-unix t))))
                     (end-pos (1+ (length (decode-coding-string
                                          (substring bytes 0 end) 'utf-8-unix t))))
                     (table (mapcar
                             (lambda (entry)
                               (propertize (car entry)
                                           'nissa-script-documentation (nth 2 entry)
                                           'nissa-script-location (nthcdr 3 entry)))
                             entries)))
                (list beg-pos end-pos table
                      :exclusive 'no
                      :exit-function
                      (lambda (candidate status)
                        (when (and (memq status '(finished exact))
                                   (assoc-string candidate entries))
                          (nissa-script-completion-insert-call (assoc-string candidate entries))))
                      :company-docsig
                      (lambda (candidate)
                        (when-let* ((entry (assoc-string candidate entries)))
                          (nissa-script-completion-summary entry)))
                      :company-doc-buffer
                      (lambda (candidate)
                        (when-let* ((entry (assoc-string candidate entries)))
                          (nissa-script-completion-doc-buffer entry)))
                      :annotation-function
                      (lambda (candidate)
                        (when-let* ((entry (assoc-string candidate entries)))
                          (concat "  "
                                  (if-let* ((signature (nissa-script-completion-signature entry)))
                                      (substring signature (length (car entry)))
                                    (nth 1 entry))
                                  (if (string-empty-p (nth 3 entry)) ""
                                    (format "  %s:%s:%s" (nth 3 entry)
                                            (nth 4 entry) (nth 5 entry)))))))))))))))

(define-minor-mode nissa-script-completion-mode
  "Complete NICE scripts through `nissa-script-completion-program'."
  :lighter " NICE"
  (if nissa-script-completion-mode
      (progn
        (add-hook 'completion-at-point-functions
                  #'nissa-script-completion-at-point nil t)
        (add-hook 'eldoc-documentation-functions #'nissa-script-completion-eldoc nil t)
        (eldoc-mode 1))
    (remove-hook 'completion-at-point-functions
                 #'nissa-script-completion-at-point t)
    (remove-hook 'eldoc-documentation-functions #'nissa-script-completion-eldoc t)
    (unless (bound-and-true-p nissa-script-check-mode)
      (nissa-script-session-stop))))

(provide 'nissa-script-completion)
;;; nissa-script-completion.el ends here
