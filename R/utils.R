emphasize <- green $ bold

# Reload method header for Verbose object to avoid too long headers
header.Verbose <- function (this, ..., char = "-", padding = 1, prefix = paste(char,
                                                                               paste(rep(" ", max(padding, 1)), collapse = ""), sep = ""),
                            level = this$defaultLevel)
{
  if (!isVisible(this, level))
    return(invisible(FALSE))
  ruler(this, char = char, length = 40)
  for (kk in seq_len(padding)) writeRaw(this, prefix, "\n")
  cat(prefix, ..., sep = "", collapse = "\n")
  for (kk in seq_len(padding)) writeRaw(this, prefix, "\n")
  ruler(this, char = char, length = 40)
}
