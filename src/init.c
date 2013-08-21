
#include <mrc.h>
#include <mrc_params.h>
#include <mrc_obj.h>

#include <assert.h>

unsigned long mrc_flags = 0;

void
mrc_set_flags(unsigned long flags)
{
  mrc_flags |= flags;
}

void
mrc_clear_flags(unsigned long flags)
{
  mrc_flags &= ~flags;
}

void
libmrc_finalize(bool assert_clean, int class_info_verbosity)
{
  int status = 0;

  struct mrc_param_select class_info_verb_options[] = {
    { .val = CLASS_INFO_VERB_NONE,   .str = "none"    },
    { .val = CLASS_INFO_VERB_DIFF,   .str = "diff"    },
    { .val = CLASS_INFO_VERB_ACTIVE, .str = "active"  },
    { .val = CLASS_INFO_VERB_FULL,   .str = "full"    },
    {},
  };
  mrc_params_get_option_select_help("class_info_verb", class_info_verb_options,
    &class_info_verbosity, "class info verbosity level, useful for leak checking");
  mrc_params_get_option_bool_help("assert_clean", &assert_clean,
    "whether to assert that things finalize cleanly");

  libmrc_params_finalize();
  status += mrc_obj_print_class_info(class_info_verbosity);

  if (assert_clean) {
    assert(status == 0);
  }
}
