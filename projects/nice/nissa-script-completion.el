;;; nissa-script-completion.el --- CAPF adapter for nASTy queries -*- lexical-binding: t; -*-

;; Backend-independent adapter for formatScriptCompletionForEmacs
;; Querying must not execute the input script

(require 'subr-x)

(defgroup nissa-script nil
  "Editing NICE scripts."
  :group 'languages)

(defcustom nissa-script-completion-program "NICE"
  "NICE executable used for completion queries."
  :type 'string
  :group 'nissa-script)

(defun nissa-script-completion-query (text byte-offset)
  "Send TEXT and BYTE-OFFSET to NICE without evaluating the script."
  (let ((errors (make-temp-file "nissa-completion-"))
        (coding-system-for-read 'utf-8-unix)
        (coding-system-for-write 'utf-8-unix))
    (unwind-protect
        (with-temp-buffer
          (insert text)
          (let ((status (call-process-region
                         (point-min) (point-max) nissa-script-completion-program
                         t (list t errors) nil
                         "--complete-at" (number-to-string byte-offset)
                         "--completion-format" "emacs")))
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
The default provider starts one NICE completion process per query, before
normal NISSA/MPI initialization, and never executes TEXT.
If using a textual response, read its Lisp data with `read', never `eval'.")

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
                                          (substring bytes 0 begin) 'utf-8 t))))
                     (end-pos (1+ (length (decode-coding-string
                                          (substring bytes 0 end) 'utf-8 t))))
                     (table (mapcar
                             (lambda (entry)
                               (propertize (car entry)
                                           'nissa-script-documentation (nth 2 entry)
                                           'nissa-script-location (nthcdr 3 entry)))
                             entries)))
                (list beg-pos end-pos table
                      :exclusive 'no
                      :annotation-function
                      (lambda (candidate)
                        (when-let* ((entry (assoc-string candidate entries)))
                          (concat "  " (nth 1 entry)
                                  (if (string-empty-p (nth 3 entry)) ""
                                    (format "  %s:%s:%s" (nth 3 entry)
                                            (nth 4 entry) (nth 5 entry)))))))))))))))

(define-minor-mode nissa-script-completion-mode
  "Complete NICE scripts through `nissa-script-completion-program'."
  :lighter " NICE"
  (if nissa-script-completion-mode
      (add-hook 'completion-at-point-functions
                #'nissa-script-completion-at-point nil t)
    (remove-hook 'completion-at-point-functions
                 #'nissa-script-completion-at-point t)))

(provide 'nissa-script-completion)
;;; nissa-script-completion.el ends here
