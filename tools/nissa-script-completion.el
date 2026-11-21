;;; nissa-script-completion.el --- CAPF adapter for nASTy queries -*- lexical-binding: t; -*-

;; Backend-independent adapter for formatScriptCompletionForEmacs
;; Querying must not execute the input script

(require 'subr-x)

(defvar-local nissa-script-completion-query-function nil
  "Function accepting source TEXT and a zero-based UTF-8 BYTE-OFFSET.
Return nil, or (BEGIN END CANDIDATES), where BEGIN and END are zero-based
UTF-8 byte offsets and each candidate is
(NAME KIND DOCUMENTATION SOURCE LINE COLUMN OFFSET).
The query must not execute TEXT.  Transport to a C++ scope provider is
deliberately supplied by the caller; this adapter starts no process.
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

(provide 'nissa-script-completion)
;;; nissa-script-completion.el ends here
