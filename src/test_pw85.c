#include <glib.h>

int main(int argc, char **argv) {
  g_test_init(&argc, &argv, NULL);

  return g_test_run();
}
